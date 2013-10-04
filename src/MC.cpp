typedef unsigned uint;
typedef double Float;
typedef long double LongFloat;

static bool verbose = false;
long double debug_num = 0;

//#define USE_RdRand

#include "util.hh"

struct SpinFlipper;
struct SpinFunc;

// random

#ifdef USE_RdRand
typedef RdRand          RandomEngine;
#else
typedef std::mt19937_64 RandomEngine;
#endif

uint64_t random_seed = 0;
RandomEngine                          random_engine;
std::uniform_int_distribution<int>    bool_dist(0, 1);
std::uniform_real_distribution<Float> uniform_dist;
std::normal_distribution<Float>       normal_dist;

// MC

class MC
{
  enum class UpdateMethod {local=0, global, smart, nMethods};
public:
  // lattice
  
  typedef uint Size;
  typedef Size Index;
  
  //
  
  static unique_ptr<MC> make(const vector<uint> &L, uint n_fields);
  
  virtual ~MC();
  
  virtual void clear_spins() =0;
  virtual void randomize_spins() =0;
  virtual void normalize_spins() =0;
  
  virtual void set_flipper(SpinFlipper *flipper) =0;
  virtual void set_potential(const SpinFunc *V) =0;
  
  virtual const SpinFunc* get_potential() const =0;
  
  void set_update_method(const string method);
  
  void sweep(const uint64_t nSweeps)
  {
    if ( _update_method == UpdateMethod::global && get_potential() )
      std::cerr << "WARNING: global updates may not be ergodic with a nonzero potential" << std::endl;
    
    FOR(sweepNum, nSweeps) {
      switch (_update_method) {
        case UpdateMethod:: local:  local_update(); break;
        case UpdateMethod::global: global_update(); break;
        case UpdateMethod:: smart: sweepNum%2 && sweepNum>0 ? local_update() : global_update(); break;
        default:
          Assert(false, int(_update_method));
      }
      
      if ( (sweepNum%100)==0 )
        normalize_spins(); // to prevent roundoff error
    }
  }
  
  Float              J;
  const Size         N; // number of spins
  const vector<uint> L_;
  const uint         n_;
  
  friend ostream& operator<<(ostream &os, const MC &mc);
  friend ostream& operator<<(ostream &os, UpdateMethod method);
  
protected:
  MC(Size N_, const vector<uint> &L_v, uint n__)
    : J(0),
      N(N_),
      L_(L_v),
      n_(n__),
     _update_method(UpdateMethod::smart),
      index_dist(0, N_-1),
     _nSweeps(0)
  { }
  
  virtual void  local_update() =0;
  virtual void   fast_update() =0;
  virtual void global_update() =0;
  
  UpdateMethod _update_method;
  
  std::uniform_int_distribution<Index> index_dist;
  
  uint64_t _nSweeps;
  vector<Float> _sum2, _sum2_2, _sum4;
};
MC::~MC() {}

ostream& operator<<(ostream &os, const MC &mc)
{
  os << "{\n"
     << "\"n\" -> " << mc.n_ << ",\n"
     << "\"L\" -> " << mc.L_ << ",\n"
     << "\"J\" -> " << mc.J << ",\n"
     << "\"update method\" -> \"" << mc._update_method << "\",\n"
  #ifndef USE_RdRand
     << "\"random seed\" -> " << random_seed << ",\n"
  #endif
     << "\"sweeps\" -> "  << mc._nSweeps << ",\n"
     << "\"s^2\" -> "     << mc._sum2    << ",\n"
     << "\"(s^2)^2\" -> " << mc._sum2_2  << ",\n"
     << "\"s^4\" -> "     << mc._sum4    << ",\n";
  return os << "}\n";
}

ostream& operator<<(ostream &os, MC::UpdateMethod method)
{
  String method_str = nullptr;
  switch (method) {
    case MC::UpdateMethod:: local: method_str =  "local"; break;
    case MC::UpdateMethod::global: method_str = "global"; break;
    case MC::UpdateMethod:: smart: method_str =  "smart"; break;
    default:
      Assert(false, int(method));
  }
  return os << method_str;
}

void MC::set_update_method(const string method_str)
{
  FOR(method_num, int(UpdateMethod::nMethods))
    if ( method_str == stringify(UpdateMethod(method_num)) ) {
      _update_method = UpdateMethod(method_num);
      return;
    }
  
  Assert(false, method_str);
}

// Spin_

template<uint n, typename _float=Float>
struct Spin_
{
  static struct Zero {} zero;

  Spin_() = default;
  Spin_(Zero) { FOR(k,n) _s[k] = 0; }
  Spin_(const Spin_ &s) = default;
  Spin_(const array<_float,n> &&s) : _s(s) {}
  
  template<typename _float_> explicit
  Spin_(const Spin_<n,_float_> &s)
  { FOR(k,n) _s[k] = s[k]; }
  
  _float& operator[](uint k)       { return _s[k]; }
  _float  operator[](uint k) const { return _s[k]; }
  
  template<typename _float_>
  void  operator+=(Spin_<n,_float_> s)       { FOR(k,n) _s[k] += s[k]; }
  void  operator-=(Spin_            s)       { FOR(k,n) _s[k] -= s[k]; }
  Spin_ operator+ (Spin_            s) const { s += _this; return s; }
  Spin_ operator- (Spin_            s) const { Spin_ ret; FOR(k,n) ret[k] =  _this[k] - s[k]; return ret; }
  Spin_ operator- (                  ) const { Spin_ ret; FOR(k,n) ret[k] = -_this[k]       ; return ret; }
  
  void  operator*=(_float x)       { FOR(k,n) _s[k] *= x; }
  void  operator/=(_float x)       { FOR(k,n) _s[k] /= x; }
  Spin_ operator* (_float x) const { Spin_ s = _this; s *= x; return s; }
  
  _float operator|(Spin_ s) const __attribute__((pure)) { _float x=0; FOR(k,n) x += _s[k]*s[k]; return x; } // dot product
  
  void normalize() { _this /= sqrt(_this|_this); }
  
  static Spin_ random()
  {
    Spin_ s;
    Float norm;
    do {
      FOR(k,n) s[k] = normal_dist(random_engine);
      norm = s|s;
    } while (norm < .01);
    s /= sqrt(norm);
    return s;
  }
  
  array<_float,n> _s;
};

template<uint n, typename _float>
ostream& operator<<(ostream &os, const Spin_<n,_float> &s)
{ return os << s._s; }

// SpinFlipper

struct SpinFlipper {
  virtual ~SpinFlipper();
  virtual void reset() =0;
};
SpinFlipper::~SpinFlipper() {}

template<uint n>
struct SpinFlipper_ : public SpinFlipper {
  typedef Spin_<n> Spin;
  virtual ~SpinFlipper_() {}
  virtual void flip(Spin &s) const =0; // must satisfy V(R.s) = V(s)
  virtual Float delta(Spin s1, const Spin s2) const // (R.s1).s2 - s1.s2 = s1.(R - 1).s2
  { Float ret = s1|s2; flip(s1); return (s1|s2) - ret; }
};

template<uint n>
struct InvertSpin_ : public SpinFlipper_<n> {
  typedef Spin_<n> Spin;
  
  virtual void  reset() final { r = Spin::random(); }
  virtual void  flip(Spin &s) const final { s += r * -2*(s|r); }
  virtual Float delta(Spin s1, Spin s2) const final { return -2*(r|s1)*(r|s2); }
  
  Spin r;
};

template<uint n>
struct SpinOperator_ : public SpinFlipper_<n> {
  typedef Spin_<n>              Spin;
  typedef function<void(Spin&)> Operator;
  
  SpinOperator_(const vector<Operator> &&operators_)
    : operators(operators_),
     _operator(nullptr),
     _dist(0, operators.size()-1)
  { }
  
  virtual void reset() final { _operator = &operators[_dist(random_engine)]; }
  virtual void flip(Spin &s) const final { (*_operator)(s); }
  
  const vector<Operator> operators;
  
private:
  const Operator *_operator;
  std::uniform_int_distribution<uint> _dist;
};

constexpr Float sqrt1_2 = sqrt(Float(.5));
SpinOperator_<1> ising_flipper({ [](Spin_<1> &s){s = -s;} });
#define OP(...) [](Spin_<COUNT_ARGS(__VA_ARGS__)> &s) { s={{__VA_ARGS__}}; }
#define r sqrt1_2
SpinOperator_<2>        permutation_flipper_2({ OP(s[1],s[0]) });
SpinOperator_<2> signed_permutation_flipper_2({ OP(-s[0],-s[1]), OP(-s[0],s[1]), OP(-s[1],-s[0]), OP(s[1],s[0]), OP(s[0],-s[1]) });
SpinOperator_<3>        permutation_flipper_3({ OP(s[2],s[1],s[0]), OP(s[1],s[0],s[2]), OP(s[0],s[2],s[1]) });
SpinOperator_<3> signed_permutation_flipper_3({ OP(-s[0],-s[1],-s[2]), OP(-s[0],-s[1],s[2]), OP(-s[0],-s[2],-s[1]), OP(-s[0],s[2],s[1]), OP(-s[0],s[1],-s[2]), OP(-s[0],s[1],s[2]), OP(-s[1],-s[0],-s[2]), OP(-s[1],-s[0],s[2]), OP(-s[2],-s[1],-s[0]), OP(-s[2],s[1],-s[0]), OP(s[2],-s[1],s[0]), OP(s[2],s[1],s[0]), OP(s[1],s[0],-s[2]), OP(s[1],s[0],s[2]), OP(s[0],-s[1],-s[2]), OP(s[0],-s[1],s[2]), OP(s[0],-s[2],-s[1]), OP(s[0],s[2],s[1]), OP(s[0],s[1],-s[2]) });
#include "SpinOperatorData.hh" // vison_triangle_flipper
#undef OP
#undef r

// SpinFunc

struct SpinFunc { virtual ~SpinFunc(); };
SpinFunc::~SpinFunc() {}

template<uint n>
struct SpinFunc_ : public SpinFunc {
  virtual ~SpinFunc_() {}
  virtual Float operator()(Spin_<n> s) const =0;
};

template<uint n, uint p=4>
struct s4_Potential_ : public SpinFunc_<n> {
  s4_Potential_(Float u_) : u(u_) {}
  virtual Float operator()(Spin_<n> s) const final { Float V=0; FOR(k,n) V += Pow<p>(s[k]); return u*V; }
  Float u;
};

unique_ptr<SpinFunc> make_s4_Potential(double u, int n, int p=4)
{
  SpinFunc *potential = nullptr;
  if        (n==2) {
    if (p==4) potential = new s4_Potential_<2,4>(u); else
    if (p==6) potential = new s4_Potential_<2,6>(u); else
    if (p==8) potential = new s4_Potential_<2,8>(u);
  } else if (n==3) {
    if (p==4) potential = new s4_Potential_<3,4>(u); else
    if (p==6) potential = new s4_Potential_<3,6>(u);
  }
  
  if ( potential )
    return unique_ptr<SpinFunc>(potential);
  else
    Assert(false, n, p);
}

struct VisonSquare_sVBS_Potential : public SpinFunc_<4>
{
  typedef Spin_<4> Spin;
  
  VisonSquare_sVBS_Potential(Float u_, Float v_) : u(u_), v(v_) {}
  
  virtual Float operator()(Spin s) const final {
    FOR(k,4) s[k] = s[k]*s[k];
    return   u*(s[0] + s[1] + s[2] + s[3])
           + v*(s[0] + s[1])*(s[2] + s[3]);
  }
  
  Float u, v;
};

struct VisonTrianglePotential : public SpinFunc_<6>
{
  typedef Spin_<6> Spin;
  
  VisonTrianglePotential(Float u_, Float v_) : u(u_), v(v_) {}
  
  virtual Float operator()(Spin s) const final {
    Spin s2;
    FOR(k,6) s2[k] = s[k]*s[k];
    return   u*(Pow<2>(s2[0]+s2[1])   +   Pow<2>(s2[2]+s2[3])   +   Pow<2>(s2[4]+s2[5]))
         + 2*v*(      (s2[0]-s2[1])*s[2]*s[3] + (s2[2]-s2[3])*s[4]*s[5] + (s2[4]-s2[5])*s[0]*s[1]);
  }
  
  Float u, v;
};

// MC_

template<uint dim, // # of spacial dimensions
         uint n>   // # of spin components
class MC_ : public MC
{
public:
  typedef Spin_<n,    Float> Spin;
  typedef Spin_<n,LongFloat> SpinSum;
  typedef array<Index,2*dim> NearestNeighbors;
  
  static_assert( std::is_pod<Spin>::value, "Spin isn't POD" );
  
  // lattice
  
  typedef array<uint,dim> Pos;
  
  Index index(const Pos p) const // todo: try blocking
  {
    Index i = p[0];
    for (uint d=1; d<dim; ++d)
      i = i*L[d] + p[d];
    return i;
  }
  
  Pos pos(Index i) const // todo: try blocking
  {
    Pos p;
    for (uint d=dim-1; d>0; --d) {
      p[d] = i%L[d];
      i /= L[d];
    }
    p[0] = i;
    return p;
  }
  
  NearestNeighbors nearestNeighbors(const Index i)
  {
    NearestNeighbors nn;
    const Pos p = pos(i);
    
    uint k = 0;
    FOR(d, dim)
    for (int dir=-1; dir<=+1; dir+=2) {
      Pos q = p;
      q[d] = (q[d] + L[d] + dir) % L[d];
      const Index j = index(q);
      
      nn[k++] = j;
    }
    
    return nn;
  }
  
  // MC
  
  virtual void clear_spins() final {
    Spin s;
    FOR(k,n) s[k] = (k==0);
    FOR(i,N) _spins[i] = s;
  }
  
  virtual void randomize_spins() final
  { FOR(i,N) _spins[i] = Spin::random(); }
  
  virtual void normalize_spins() final
  { FOR(i,N) _spins[i].normalize(); }
  
  virtual void local_update() final __attribute__((hot))
  {
    for (Size count=0; count<N; ++count) { // warning: do NOT count how many spins are updated
      const Index i = index_dist(random_engine);
      
      const Spin s1 = _spins[i];
      const Spin s2 = n>1 ? Spin::random() : -s1;
      const Spin ds = s2 - s1;
      Float delta_E = 0;
      if (_V)
        delta_E += V(s2) - V(s1);
      for(Index nn : nearestNeighbors(i))
        delta_E += -J*(ds|_spins[nn]);
      
      if ( delta_E <= 0 || uniform_dist(random_engine) < exp(-delta_E) )
        _spins[i] = s2;
    }
    
    measure();
    
    ++_nSweeps;
    _update_sublattice = 2;
  }
  
  virtual void fast_update() final __attribute__((hot))
  {
    FOR(d, dim)
      Assert( L[d]%2 == 0, L[d] );
    
    if ( _update_sublattice == 2 )
      _update_sublattice = bool_dist(random_engine);
    else
      _update_sublattice = !_update_sublattice;
    
    for (Index i0=0; i0<N; i0+=N/L[dim-1]) {
      const Pos p0 = pos(i0);
      Index i = _update_sublattice;
      FOR(d, dim-1)
        i += p0[d];
      
      for (i = i%2; i<i0+L[dim-1]; i+=2) {
        //const Spin s1 = _spins[i];
        Spin sn(Spin::zero);
        
        for(Index nn : nearestNeighbors(i))
          sn += _spins[nn];
        
        //TODO
      }
    }
    
    ++_nSweeps;
  }
  
  virtual void global_update() final __attribute__((hot))
  {
    if ( _check_flipper )
    {
      FOR(i, 100) {
        _flipper->reset();
        Spin s = Spin::random();
        const Float V1 = V(s);
        _flipper->flip(s);
        const Float V2 = V(s);
        Assert( abs(V1-V2) < sqrt(std::numeric_limits<Float>::epsilon()), s, V1, V2 );
      }
      _check_flipper = false;
    }
    
    _flipper->reset();
    _cluster.assign(N, false);
    
    Size nClusters = 0;
    FOR(i0, N)
    if (!_cluster[i0]) {
      const bool flipCluster = bool_dist(random_engine);
      Index *newIndex = _newIndexStack.get();
      
      *newIndex = i0;
      _cluster[i0] = true;
      //_clusterSums[nClusters] = SpinSum(_spins[i0]);
      if ( flipCluster )
        _flipper->flip(_spins[i0]);
      
      do {
        const Index j = *(newIndex--);
        
        for(Index i : nearestNeighbors(j))
          if ( !_cluster[i] ) {
            const Float delta_E = -J*(1-2*flipCluster) * _flipper->delta(_spins[i], _spins[j]);
            //debug_num += exp(uniform_dist(random_engine));
            if ( delta_E > 0 && uniform_dist(random_engine) > exp(-delta_E) ) {
              *(++newIndex) = i;
              _cluster[i] = true;
              //_clusterSums[nClusters] += _spins[i];
              if ( flipCluster )
                _flipper->flip(_spins[i]);
            }
          }
      }
      while ( newIndex+1 != _newIndexStack.get() );
      ++nClusters;
    }
    /*
    SpinSum X(SpinSum::zero); // spin sum part that doesn't get flipped
    SpinSum Y(SpinSum::zero), // spin sum part that gets flipped
            long_r(r);
    LongFloat Y2=0, Y4=0;
    
    FOR(c, nClusters) {
      const SpinSum spin = _clusterSums[c];
      const LongFloat y = spin | long_r;
      X  += spin - long_r*y;
      Y2 += Pow<2>(y);
      Y4 += Pow<4>(y);
    }
    const LongFloat
      X2  = X|X,
      Xr2 = Pow<2>(X|long_r),
      Z2  = Y2,
      Z4  = 3*Y2 - 2*Y4;
    
    //_sum2.push_back( X2 + Z2 );
    //_sum4.push_back( Pow<2>(X2) + (2*X2 + Xr2)*Z2 + Z4 );
    */
    measure(); // TODO
    
    ++_nSweeps;
    _update_sublattice = 2;
  }
  
  virtual void set_flipper(SpinFlipper *flipper) final
  {
    if (flipper) {
      _flipper = dynamic_cast<SpinFlipper_<n>*>(flipper);
      Assert(_flipper, n);
    }
    else
      _flipper = n>1 ? &_invert_spin : dynamic_cast<SpinFlipper_<n>*>(n==1 ? &ising_flipper : nullptr);
    
    if (_V)
      _check_flipper = true;
  }
  
  virtual void set_potential(const SpinFunc *new_V) final
  {
    _V = dynamic_cast<const SpinFunc_<n>*>(new_V);
    
    if (new_V) {
      _check_flipper = true;
      Assert(_V, n);
    }
  }
  
  virtual const SpinFunc* get_potential() const final { return _V; }
  
  Float V(Spin s) const { return (*_V)(s); }
  
protected:
  void measure()
  {
    SpinSum sum(SpinSum::zero);
    FOR(i, N)
      sum += _spins[i];
    const LongFloat sum2 = sum|sum;
    LongFloat sum4  = 0;
    FOR(k,n)  sum4 += Pow<4>(sum[k]);
    
    _sum2  .push_back(sum2);
    _sum2_2.push_back(sum2*sum2);
    _sum4  .push_back(sum4);
  }
  
  MC_(const array<uint,dim> &L__, const vector<uint> &L_v, Size N_)
    : MC(N_, L_v, n),
      L(L__),
     _spins(new Spin[N_]),
     _newIndexStack(new Index[N_]),
     _clusterSums(new SpinSum[N_]),
     _flipper(nullptr),
     _V(nullptr),
     _update_sublattice(2),
     _check_flipper(false)
  { set_flipper(nullptr); }
  
  friend class MC;
  
  virtual ~MC_() {}
  
  const array<uint,dim> L; // lengths
  
private:
  const unique_ptr<Spin[]>     _spins;
  vector<bool>                 _cluster;
  const unique_ptr<Index[]>    _newIndexStack;
  const unique_ptr<SpinSum[]>  _clusterSums;
  
  InvertSpin_<n>               _invert_spin;
  SpinFlipper_<n>             *_flipper;
  const SpinFunc_<n>          *_V;
  
  uint _update_sublattice; // 0: A sublattice, 1: B sublattice, 2: random
  bool _check_flipper;
};

unique_ptr<MC> MC::make(const vector<uint> &L, uint n_fields)
{
  MC::Size N = 1;
  const uint dim = L.size();
  FOR(d,dim) {
    const MC::Size old_N = N;
    N *= L[d];
    Assert( L[d] && N/L[d] == old_N, old_N, L[d], N ); // check overflow
  }
  
  MC *mc = nullptr;
  #define ELSE_TRY_MC(dim_, n_) \
  else if (dim == dim_ && n_fields == n_) \
  { \
    array<uint,dim_> L_; \
    FOR(d,dim) L_[d] = L[d]; \
    mc = new MC_<dim_,n_>(L_, L, N); \
  }
  
  if (false) {}
//ELSE_TRY_MC(1,1)
//ELSE_TRY_MC(1,2)
//ELSE_TRY_MC(2,1)
//ELSE_TRY_MC(2,2)
//ELSE_TRY_MC(2,3)
//ELSE_TRY_MC(3,1)
  ELSE_TRY_MC(3,2)
  ELSE_TRY_MC(3,3)
  ELSE_TRY_MC(3,4)
  ELSE_TRY_MC(3,6)
//ELSE_TRY_MC(4,1)
//ELSE_TRY_MC(4,2)
//ELSE_TRY_MC(4,3)
  else
    Assert(false, L, N);
  
  #undef ELSE_TRY_MC
  
  return unique_ptr<MC>(mc);
}

static void single_handler(int sig)
{
  printf("Error: signal %d:\n", sig);
  
  print_backtrace();
  
  exit(1);
}

int main(const int argc, char *argv[])
{
  using std::cout;
  using std::cerr;
  using std::endl;
  
  // http://www.gnu.org/s/hello/manual/libc/Signal-Handling.html
  signal( SIGSEGV, single_handler );
  //signal( SIGABRT, single_handler );
  
  // read program options
  namespace po = boost::program_options;
  
  po::options_description generic_options("generic options");
  generic_options.add_options()
      ("help,h",                               "print help message")
      ("verbose,v",                            "print verbose output");
  
  uint n = 0;
  vector<double> J;
  string potential_str;
  po::options_description system_options("physics options");
  system_options.add_options()
      ("L",         po::value<vector<uint>>()->multitoken(),     "lengths")
      ("n",         po::value<uint>(&n),                         "for an O(n) model")
      ("J",         po::value<vector<double>>(&J)->multitoken(), "spin coupling and other potential coefficients")
      ("potential", po::value<string>(&potential_str),           "potential term: s^n where n=4,6,8; vison hexagon|square|triangle c|s-VBS");
  
  po::options_description simulation_options("simulation options");
  simulation_options.add_options()
      ("sweeps",                 po::value<uint64_t>()->default_value(1),      "# of MC sweeps")
      ("file",                   po::value<string>(),                          "save file")
      ("update-method",          po::value<string>()->default_value("smart"),  "update type: local, global, or smart");
  
  po::options_description cmdline_options;
  cmdline_options.add(generic_options)
                 .add(system_options)
                 .add(simulation_options);
  
  po::variables_map vm;
  store(po::parse_command_line(argc, argv, cmdline_options), vm);
  notify(vm);
  
  verbose = vm.count("verbose");
  
  /// generic options
  if ( vm.count("help") ) {
    cout << cmdline_options << endl;
    return 0;
  }
  
  if ( !vm.count("L") ) {
    cerr << "L is required" << endl;
    return 1;
  }
  
  unique_ptr<SpinFunc>   potential;
  SpinFlipper          *spin_flipper = nullptr;
  if ( !potential_str.empty() )
  {
    auto check_n = [&n](const uint new_n) {
      if (n==0) n=new_n;
      else      Assert(n==new_n, n, new_n); };
    
    if        ( potential_str == "vison hexagon c-VBS" ) {
      check_n(2);
      potential_str = "s^6";
    } else if ( potential_str == "vison hexagon s-VBS" ) {
      check_n(3);
      potential_str = "s^4";
    } else if ( potential_str == "vison square c-VBS" ) {
      check_n(2);
      potential_str = "s^8";
    }
    
    if ( potential_str.substr(0,potential_str.size()-1) == "s^" ) {
      const uint pow = from_string<uint>(potential_str.substr(2));
      Assert(J.size() == 1+1, J);
      potential = make_s4_Potential(J[1], n, pow);
      if (n==2) spin_flipper = &signed_permutation_flipper_2;
      if (n==3) spin_flipper = &signed_permutation_flipper_3;
      else Assert(false, n);
    } else if ( potential_str == "vison square s-VBS" ) {
      check_n(4);
      Assert(J.size() == 1+2, J);
      potential.reset(new VisonSquare_sVBS_Potential(J[1], J[2]));
      Assert(false, potential_str);
    } else if ( potential_str == "vison triangle c-VBS" ) {
      check_n(6);
      Assert(J.size() == 1+2, J);
      potential.reset(new VisonTrianglePotential(J[1], J[2]));
      spin_flipper = &vison_triangle_flipper;
    } else
      Assert(false, potential_str);
  }
  else
    Assert(J.size() == 1, J);
  
  Assert(n, n);
  unique_ptr<MC> mc = MC::make(vm["L"].as<vector<uint>>(), n);
  
  #ifndef USE_RdRand
  random_seed = std::chrono::system_clock::now().time_since_epoch().count() * getpid();
  random_engine.seed(random_seed);
  #endif
  
  mc->clear_spins();
  mc->J = J[0];
  mc->set_update_method( vm["update-method"].as<string>() );
  mc->set_potential(potential.get());
  mc->set_flipper(spin_flipper);
  
  mc->sweep( vm["sweeps"].as<uint64_t>() );
  
  if ( vm.count("file") ) {
    std::ofstream file( vm["file"].as<string>() );
    file << *mc;
  }
  
  if ( debug_num )
    cout << debug_num << endl;
  
  return 0;
}
