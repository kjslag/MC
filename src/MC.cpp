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

// Value

struct Value
{
  Value() : _n_samples(0), _means_bin_num(0) { _means.reserve(max_means); }
  
  void operator<<(LongFloat x)
  {
    Assert( _means_bin_num == calc_means_bin_num(),
            max_means, _n_samples, _means_bin_num, calc_means_bin_num() );
    
    uint k = 0; // bin number
    next_bin: {
      if ( k == _means_bin_num ) {
        if ( _means.size() < max_means )
          _means.push_back(x/exp2i(_means_bin_num));
        else {
          FOR(i, max_means/2)
            _means[i] = (_means[2*i] + _means[2*i+1])/2;
          _means.resize(max_means/2);
          Assert( !full(k), k, _n_samples ); // else there would be double counting
          ++_means_bin_num;
           // x will be added later
        }
      }
      
      if ( k < _binned_vec.size() ) {
        Binned &binned = _binned_vec[k];
        binned.x2 += x*x;
        if ( !full(k) )
          binned.x = x;
        else {
          x += binned.x;
          binned.x = 0;
          ++k;
          goto next_bin;
        }
      } else
        _binned_vec .push_back(x);
    }
    
    ++_n_samples;
  }
  
  uint default_error_bin_num() const
  { return max(int(_binned_vec.size()) - 6, 0); }
  
  LongFloat mean() const
  {
    // TODO: calc from _binned_vec
    const uint n_means = _means.size();
    LongFloat mean_ = 0;
    if ( n_means >= 2 ) {
      for (uint i=1; i<n_means; ++i)
        mean_ += _means[i];
      mean_ /= n_means - 1;
    } else
      mean_ = NAN;
    
    if ( _n_samples>1 ) {
      const uint bin_num = default_error_bin_num();
      LongFloat x = mean(bin_num);
      Assert( is_small((x-mean_)/(x+mean_)), mean_, x );
    }
    
    return mean_;
  }
  
  LongFloat mean(uint bin_num) const
  {
    LongFloat mean_ = 0;
    if ( _n_samples > 1 ) {
      const uint n_bins = _binned_vec.size();
      FOR(k, n_bins)
        mean_ += _binned_vec[k].x;
      mean_ -= _binned_vec[bin_num].first_x;
      mean_ /= error_bin_count(bin_num);
    } else
      mean_ = NAN;
    return mean_;
  }
  
  const vector<LongFloat>& means() const { return _means; }
  
  LongFloat error() const { return _n_samples>1 ? error(default_error_bin_num()) : INFINITY; }
  
  LongFloat error(uint bin_num) const
  {
    const LongFloat mean_ = mean(bin_num);
    const LongFloat x2 = _binned_vec[bin_num].x2 - _binned_vec[bin_num].first_x2;
    const uint64_t n_samples = error_bin_count(bin_num);
    return sqrt(max(LongFloat(0), x2 - n_samples*bin_size(bin_num)*mean_*mean_)) / n_samples;
  }
  
  vector<LongFloat> binned_errors() const
  {
    vector<LongFloat> errors_;
    const int n_bins = _binned_vec.size();
    FOR(bin_num, n_bins-1)
      errors_.push_back( error(bin_num) );
    return errors_;
  }
  
protected:
  struct Binned {
    Binned(LongFloat y) : x(y), x2(y*y), first_x(y), first_x2(y*y) {}
    
    LongFloat x;  // a partial sum
    LongFloat x2; // binned x^2
    const LongFloat first_x, first_x2;
  };
  
  friend ostream& operator<<(ostream &os, const Binned &binned);
  friend ostream& operator<<(ostream &os, const Value  &value);
  
  bool full(uint bin_num) const
  { return _n_samples & exp2i(bin_num); }
  
  // number of samples per bin
  static uint64_t bin_size(uint bin_num)
  { return exp2i(bin_num); }
  
  // number of samples that went into _binned_vec[bin_num].x2
  uint64_t error_bin_count(uint bin_num) const
  { return (_n_samples & ~(bin_size(bin_num)-1)) - bin_size(bin_num); }
  
  uint calc_means_bin_num() const {
    uint n = 0;
    while ( exp2i(n)*max_means < _n_samples )
      ++n;
    return n;
  }
  
private:
  uint64_t          _n_samples;
  vector<Binned>    _binned_vec; // indexed by bin_num
  vector<LongFloat> _means;
  uint              _means_bin_num;
  
  static const uint max_means;
};
const uint Value::max_means = 1u<<4; // should be a power of 2 // TODO change to 1<<8

ostream& operator<<(ostream &os, const Value::Binned &binned)
{ return os << "binned[" << binned.x << ", " << binned.x2 << ", " << binned.first_x << ", " << binned.first_x2 << "]"; }

ostream& operator<<(ostream &os, const Value &value)
{
  return os << "value["                  << value.mean()
            << ", "                      << value.error()
            << ", \"means\" -> "         << value.means()
            << ", \"binned errors\" -> " << value.binned_errors()
            << ", \"# samples\" -> "     << value._n_samples
            << ", \"_binned_vec\" -> "   << value._binned_vec
            << "]";
}

// MC

class MC
{
  enum class UpdateMethod {local=0, global, nMethods};
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
  void         set_potential(const SpinFunc *V, string potential_name) {
    _potential_name = potential_name;
    set_potential_func(V);
  }
  
  virtual const SpinFunc* potential() const =0;
  
  void set_update_method(const string method);
  
  void sweep(const uint64_t nSweeps)
  {
    FOR(sweepNum, nSweeps) {
      switch (_update_method) {
        case UpdateMethod:: local:
          local_update();
          break;
        case UpdateMethod::global:
          if (n_>1) {
            global_update(false);
            local_update (false);
          }
          global_update(true);
          break;
        default:
          Assert(false, _update_method);
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
     _update_method(UpdateMethod::global),
      index_dist(0, N_-1),
     _nSweeps(0)
  { }
  
  virtual void set_potential_func(const SpinFunc *V) =0;
  
  virtual void  local_update(bool measure=true) =0;
  virtual void global_update(bool measure) =0;
  
  UpdateMethod _update_method;
  
  std::uniform_int_distribution<Index> index_dist;
  
  string   _potential_name;
  uint64_t _nSweeps;
  vector<Float> _sum2, _sum2_2, _sum4;
};
MC::~MC() {}

ostream& operator<<(ostream &os, MC::UpdateMethod method)
{
  String method_str = nullptr;
  switch (method) {
    case MC::UpdateMethod:: local: method_str =  "local"; break;
    case MC::UpdateMethod::global: method_str = "global"; break;
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
  void  operator/=(_float x)       { _this *= 1/x; }
  Spin_ operator* (_float x) const { Spin_ s = _this; s *= x; return s; }
  
  _float operator|(Spin_ s) const { _float x=0; FOR(k,n) x += _s[k]*s[k]; return x; } // dot product
  
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
  
  void flip(Spin &s, bool measure) const { if (measure) s=-s; else flip(s); }
  Spin flipped(Spin s, bool measure) const { flip(s,measure); return s; };
  
protected:
  virtual void flip(Spin &s) const =0; // must satisfy V(R.s) = V(s), and is assumed to be linear
};

template<uint n>
struct InvertSpin_ : public SpinFlipper_<n> {
  typedef Spin_<n> Spin;
  
  virtual void  reset() final { r = Spin::random(); }
  virtual void  flip(Spin &s) const final { s += r*(-2*(s|r)); }
  
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

#include "SpinOperatorData.hh"

// SpinFunc

struct SpinFunc {
  virtual ~SpinFunc();
  virtual vector<Float> coefficients() const =0;
};
SpinFunc::~SpinFunc() {}

template<uint n>
struct SpinFunc_ : public SpinFunc {
  virtual ~SpinFunc_() {}
  virtual Float operator()(Spin_<n> s) const =0;
};

template<uint n>
struct s4_Potential_ : public SpinFunc_<n> {
  s4_Potential_(Float u_) : u(u_) {}
  virtual Float operator()(Spin_<n> s) const final { Float V=0; FOR(k,n) V += Pow<4>(s[k]); return u*V; }
  virtual vector<Float> coefficients() const { return {u}; }
  const Float u;
};

unique_ptr<SpinFunc> make_s4_Potential(Float u, uint n)
{
  SpinFunc *potential = nullptr;
  if      (n==2) potential = new s4_Potential_<2>(u);
  else if (n==3) potential = new s4_Potential_<3>(u);
  
  Assert(potential, n);
  return unique_ptr<SpinFunc>(potential);
}

template<uint m>
struct cos_Potential_ : public SpinFunc_<2> {
  cos_Potential_(Float u_) : u(u_) {}
  virtual Float operator()(Spin_<2> s) const final { complex<Float> z(s[0], s[1]); return 2*real(Pow<m>(z)); }
  virtual vector<Float> coefficients() const { return {u}; }
  const Float u;
};

unique_ptr<SpinFunc> make_cos_Potential(Float u, uint m)
{
  SpinFunc *potential = nullptr;
  if      (m==2) potential = new cos_Potential_<2>(u);
  else if (m==3) potential = new cos_Potential_<3>(u);
  else if (m==4) potential = new cos_Potential_<4>(u);
  else if (m==6) potential = new cos_Potential_<6>(u);
  else if (m==8) potential = new cos_Potential_<8>(u);

  Assert(potential, m);
  return unique_ptr<SpinFunc>(potential);
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
  
  virtual vector<Float> coefficients() const { return {u,v}; }
  
  const Float u, v;
};

struct VisonTrianglePotential : public SpinFunc_<6>
{
  typedef Spin_<6> Spin;
  
  VisonTrianglePotential(Float u_, Float v_) : u(u_), v(v_) {}
  
  virtual Float operator()(Spin s) const final {
    Spin s2;
    FOR(k,6) s2[k] = s[k]*s[k];
    return   u*(sq(s2[0]+s2[1])     +     sq(s2[2]+s2[3])     +     sq(s2[4]+s2[5]))
         + 2*v*(  (s2[0]-s2[1])*s[2]*s[3] + (s2[2]-s2[3])*s[4]*s[5] + (s2[4]-s2[5])*s[0]*s[1]);
  }
  
  virtual vector<Float> coefficients() const { return {u,v}; }
  
  const Float u, v;
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
  
  virtual void local_update(bool measure_=true) final __attribute__((hot))
  {
    for (Size count=0; count<N; ++count) { // warning: do NOT count how many spins are updated
      const Index i = index_dist(random_engine);
      
      const Spin s1 = _spins[i];
      const Spin s2 = n>1 ? Spin::random() : -s1;
      Float delta_E = 0;
      
      if (_V)
        delta_E += V(s2) - V(s1);
      
      Spin sum_nn(Spin::zero);
      for(Index nn : nearestNeighbors(i))
        sum_nn += _spins[nn];
      delta_E += -J*((s2-s1)|sum_nn);
      
      if ( delta_E <= 0 || uniform_dist(random_engine) < exp(-delta_E) )
        _spins[i] = s2;
    }
    
    if (measure_)
      measure();
    
    ++_nSweeps;
  }
  
  virtual void global_update(bool measure_) final __attribute__((hot))
  {
    if ( _check_flipper )
    {
      FOR(i, 100) {
        _flipper->reset();
        const Spin s = Spin::random();
        const Float delta_V = V(s) - V(_flipper->flipped(s,i<3));
        Assert( is_small(delta_V), s, delta_V );
      }
      _check_flipper = false;
    }
    
    if ( !measure_ )
      _flipper->reset();
    _cluster.assign(N, false);
    
    Size nClusters = 0;
    FOR(i0, N)
    if (!_cluster[i0]) {
      const bool flipCluster = bool_dist(random_engine);
      Index *newIndex = _newIndexStack.get();
      
      *newIndex = i0;
      _cluster[i0] = true;
      if ( measure_ )
        _clusterSums[nClusters] = SpinSum(_spins[i0]);
      if ( flipCluster )
        _flipper->flip(_spins[i0], measure_);
      
      do {
        const Index j = *(newIndex--);
        
        for(Index i : nearestNeighbors(j))
          if ( !_cluster[i] ) {
            const Spin s = _spins[i];
            const Spin flipped_s = _flipper->flipped(s, measure_);
            const Float delta_E = -J*(1-2*flipCluster) * ((flipped_s-s)|_spins[j]);
            if ( delta_E > 0 && uniform_dist(random_engine) > exp(-delta_E) ) {
              *(++newIndex) = i;
              _cluster[i] = true;
              if ( measure_ )
                _clusterSums[nClusters] += s;
              if ( flipCluster )
                _spins[i] = flipped_s;
            }
          }
      }
      while ( newIndex+1 != _newIndexStack.get() );
      ++nClusters;
    }
    
    if ( measure_ ) {
      if (true) {
        LongFloat Y2=0, Y4=0;
        FOR(c, nClusters) {
          const SpinSum y  = _clusterSums[c];
          const Float   y2 = y|y;
          Y2 += y2;
          Y4 += y2*y2;
        }
        const LongFloat
          Z2 = Y2,
          Z4 = 3*sq(Y2) - 2*Y4;
        
        _sum2  .push_back(( Z2 )/sq(LongFloat(N)));
        _sum2_2.push_back(( Z4 )/Pow<4>(LongFloat(N)));
        // TODO sum4
      } else
        measure();
      
      ++_nSweeps;
    }
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
  
  virtual const SpinFunc* potential() const final { return _V; }
  
  Float V(Spin s) const { return (*_V)(s); }
  
protected:
  void measure()
  {
    SpinSum sum(SpinSum::zero);
    FOR(i, N)
      sum += _spins[i];
    sum /= N;
    const LongFloat sum2 = sum|sum;
    LongFloat sum4  = 0;
    FOR(k,n)  sum4 += Pow<4>(sum[k]);
    
    _sum2  .push_back(sum2);
    _sum2_2.push_back(sum2*sum2);
    if (n > 1)
      _sum4.push_back(sum4);
  }
  
  MC_(const array<uint,dim> &L__, const vector<uint> &L_v, Size N_)
    : MC(N_, L_v, n),
      L(L__),
     _spins(new Spin[N_]),
     _newIndexStack(new Index[N_]),
     _clusterSums(new SpinSum[N_]),
     _flipper(nullptr),
     _V(nullptr),
     _check_flipper(false)
  { set_flipper(nullptr); }
  
  virtual void set_potential_func(const SpinFunc *new_V) final
  {
    _V = dynamic_cast<const SpinFunc_<n>*>(new_V);
    
    if (new_V) {
      _check_flipper = true;
      Assert(_V, n);
    }
  }
  
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
  ELSE_TRY_MC(2,1)
  ELSE_TRY_MC(2,2)
  ELSE_TRY_MC(2,3)
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

ostream& operator<<(ostream &os, const MC &mc)
{
  const SpinFunc *const potential = mc.potential();
  vector<Float> J{mc.J};
  if (potential) {
    const vector<Float> coefficients = potential->coefficients();
    J.insert(J.end(), coefficients.begin(), coefficients.end());
  }
  
  os << "{\n"
        "\"L\" -> " << mc.L_ << ",\n"
        "\"n\" -> " << mc.n_ << ",\n"
        "\"J\" -> " <<    J  << ",\n"
        "\"potential\" -> \""     << mc._potential_name << "\",\n"
        "\"update method\" -> \"" << mc._update_method << "\",\n"
         "\"sweeps\" -> "         << mc._nSweeps << ",\n"
            "\"s^2\" -> value@"   << mc._sum2    << ",\n"
        "\"(s^2)^2\" -> value@"   << mc._sum2_2  << ",\n"
            "\"s^4\" -> value@"   << mc._sum4    << "\n";
  return os << "}\n";
}

int main(const int argc, char *argv[])
{
  using std::cout;
  using std::cerr;
  using std::endl;
  
  Value v;
  FOR(i,200) {
    cout << v << endl;
    v << 1;
  }
  
  return 0;
  
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
  vector<Float> J;
  string potential_name;
  po::options_description system_options("physics options");
  system_options.add_options()
      ("L",         po::value<vector<uint>>()->multitoken(),     "lengths")
      ("n",         po::value<uint>(&n),                         "for an O(n) model")
      ("J",         po::value<vector<Float>>(&J)->multitoken(),  "spin coupling and other potential coefficients")
      ("potential", po::value<string>(&potential_name),          "potential term: s^4; cos#; vison hexagon|square|triangle c|s-VBS");
  
  po::options_description simulation_options("simulation options");
  simulation_options.add_options()
      ("sweeps",                 po::value<uint64_t>()->default_value(1),      "# of MC sweeps")
      ("file",                   po::value<string>(),                          "save file")
      ("update-method",          po::value<string>()->default_value("global"), "update type: local or global");
  
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
  if ( !potential_name.empty() )
  {
    auto check_n = [&n](const uint new_n) {
      if (n==0) n=new_n;
      else      Assert(n==new_n, n, new_n); };
    
    string V_str = potential_name;
    
    if        ( V_str == "vison hexagon c-VBS" ) {
      V_str = "cos6";
    } else if ( V_str == "vison hexagon s-VBS" ) {
      check_n(3);
      V_str = "s^4";
    } else if ( V_str == "vison square c-VBS" ) {
      V_str = "cos8";
    }
    
    if ( V_str == "s^4" ) {
      Assert(J.size() == 1+1, J);
      potential = make_s4_Potential(J[1], n);
      if      (n==2) Assert(false, n); //spin_flipper = &signed_permutation_flipper_2;
      else if (n==3) spin_flipper = &signed_permutation_flipper_3;
      else if (n==4) spin_flipper = &signed_permutation_flipper_4;
      else           Assert(false, n);
    } else if ( V_str.substr(0,3) == "cos" ) {
      check_n(2);
      Assert(J.size() == 1+1, J);
      const uint m = from_string<uint>(V_str.substr(3));
      potential = make_s4_Potential(J[1], m);
      if      (m==2) spin_flipper = &cos2_flipper;
      else if (m==4) spin_flipper = &cos4_flipper;
      else if (m==6) spin_flipper = &cos6_flipper;
      else if (m==8) spin_flipper = &cos8_flipper;
      else           Assert(false, m);
    } else if ( V_str == "vison square s-VBS" ) {
      check_n(4);
      Assert(J.size() == 1+2, J);
      potential.reset(new VisonSquare_sVBS_Potential(J[1], J[2]));
      Assert(false, V_str);
    } else if ( V_str == "vison triangle c-VBS" ) {
      check_n(6);
      Assert(J.size() == 1+2, J);
      potential.reset(new VisonTrianglePotential(J[1], J[2]));
      spin_flipper = &vison_triangle_flipper;
    } else
      Assert(false, V_str);
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
  Assert(J.size(), J);
  mc->J = J[0];
  mc->set_update_method( vm["update-method"].as<string>() );
  mc->set_potential(potential.get(), potential_name);
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
