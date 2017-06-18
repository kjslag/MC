// developer: Kevin Slagle (kslagle@physics.ucsb.edu)
// date: 2013

typedef unsigned uint;
typedef double Float;
typedef long double LongFloat;

static bool __verbose = false;
long double debug_num = 0;

//#define USE_RdRand

#include "util.hh"

typedef IO_float<Float> IOFloat;

const Float pi = boost::math::constants::pi<Float>();

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
  Value() : Value(0,0,0) {}
  
  // min and max histo values
  // the add() function must be used instead of operator<<
  // be careful not to add improved (reduced variance) estimators to the histo
  Value(LongFloat min, LongFloat max, size_t bins=1u<<8)
   : _n_samples(0), _means_k(0),
     _histo(bins),
     _histo_min(min), _histo_max(max),
     _min_sample(+infinity<LongFloat>()), _max_sample(-infinity<LongFloat>())
  { _means.reserve(max_n_means); }
  
  void add(LongFloat x, LongFloat xh)
  {
    add_(x);
    
    if ( xh < _min_sample ) _min_sample = xh;
    if ( xh > _max_sample ) _max_sample = xh;
    
    int n = _histo.size();
    Assert(n);
    LongFloat bin = n*(xh-_histo_min)/(_histo_max-_histo_min);
    ++_histo[max(0,min<int>( n-1, lround(bin-.5) ))];
  }
  
  void add(LongFloat x) { add(x, x); }
  
  void operator<<(LongFloat x)
  { Assert(!_histo.size()); add_(x); }
  
  uint64_t n_samples() const { return _n_samples; }
  
  LongFloat mean() const
  {
    LongFloat mean_ = 0;
    if ( _n_samples ) {
      const uint n = _sums.size();
      FOR(k, n)
        mean_ += _sums[k].x;
      mean_ /= _n_samples;
    } else
      mean_ = NAN;
    
    return mean_;
  }
  
  const vector<LongFloat>& means() const { return _means; }
  
  uint default_bin_size_exp() const { return _means_k; }
  
  LongFloat error() const { return error(default_bin_size_exp()); }
  
  LongFloat error(uint k) const
  {
    if ( k+1 >= _sums.size() )
      return INFINITY;
    
    const LongFloat mean2         = sq(mean());
    const uint64_t  n_x2_samples_ = n_x2_samples(k);
    const LongFloat var           = _sums[k].x2 - n_x2_samples_*(exp2i(k)*mean2);
    return sqrt(max(var,LongFloat(0)) + 2*mean2) / n_x2_samples_;
    // + 2*mean2 in case the algorithm isn't very ergodic
  }
  
  vector<LongFloat> binned_errors() const
  {
    vector<LongFloat> errors_;
    const int n_sums = _sums.size();
    FOR(k, n_sums-1)
      errors_.push_back( error(k) );
    return errors_;
  }
  
  static bool brief;
  
protected:
  struct Sum {
    Sum(LongFloat x_) : x(x_), x2(x_*x_) {}
    
    LongFloat x;
    LongFloat x2;
  };
  
  void add_(LongFloat x)
  {
    Assert_( _means_k == calc_means_k(), _means_k );
    
    uint k = 0;
    next_bin:
    {
      if ( k == _means_k )
      {
        if ( _means.size() < max_n_means )
          _means.push_back(x/exp2i(_means_k));
        else {
          FOR(i, max_n_means/2)
            _means[i] = (_means[2*i] + _means[2*i+1])/2;
          _means.resize(max_n_means/2);
          Assert_( !full(k), k, n_samples() ); // else there would be double counting
          ++_means_k;
           // x will be added later
        }
      }
      
      if ( k < _sums.size() )
      {
        Sum &sum = _sums[k];
        sum.x2 += x*x;
        if ( !full(k) )
          sum.x = x;
        else {
          x += sum.x;
          sum.x = 0;
          ++k;
          goto next_bin;
        }
      } else
        _sums.push_back(x);
    }
    
    ++_n_samples;
  }
  
  bool full(uint k) const { return _n_samples & exp2i(k); }
  
  friend ostream& operator<<(ostream &os, const Value &value);
  
  // number of samples that went into _sums[k].x2
  uint64_t n_x2_samples(uint k) const
  { return n_samples() & ~(exp2i(k)-1); }
  
  uint calc_means_k() const {
    uint n = 0;
    while ( max_n_means < n_samples()/exp2i(n) )
      ++n;
    return n;
  }
  
private:
  uint64_t           _n_samples;
  vector<Sum>        _sums; // indexed by k, the bin size exponent
  vector<LongFloat>  _means;
  uint               _means_k;
  
  vector<uint64_t>   _histo;
  const LongFloat    _histo_min,  _histo_max;
  LongFloat          _min_sample, _max_sample;
  
  static const uint max_n_means;
};
bool       Value::brief       = false;
const uint Value::max_n_means = 1u<<6;

ostream& operator<<(ostream &os, const Value &value)
{
  os << "value["                   << value.mean()
     << ", "                       << value.error()
     << ", \"samples\" -> "        << value.n_samples();
  if ( !Value::brief ) {
  os << ",\n\"means\" -> "         << value.means()
     << ",\n\"binned errors\" -> " << value.binned_errors();
  if ( value._histo.size() ) {
  os << ",\n\"histo\" -> "         << value._histo
     << ",\n\"histo min\" -> "     << value._histo_min
     << ", \"histo max\" -> "      << value._histo_max
     << ", \"min sample\" -> "     << value._min_sample
     << ", \"max sample\" -> "     << value._max_sample;
  }}
  return os << "]";
}

// MC

struct SpinFlipper;
struct SpinFunc;

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
  void         set_potential(SpinFunc *V, string potential_name) {
    _potential_name = potential_name;
    set_potential_func(V);
  }
  
  virtual const SpinFunc* potential() const =0;
  
  void set_layer_dims(uint layer_dims) {
    _layer_dims = layer_dims;
    for (uint d=0; d<layer_dims; ++d)
      J_anisotropy[d] = 0;
  }
  void set_update_method(const string method);
  void set_thermalization(uint64_t thermalization) { _thermalization = thermalization; }
  void set_annealing(uint64_t annealing) { _annealing = annealing; }
  void set_sweeps(uint64_t nSweeps) { _n_sweeps = nSweeps; }
  
  void sweep()
  {
    Assert(_n_sweeps);
    
    #define TIME(t, f) do{ \
      using namespace std::chrono; \
      const auto t0 = steady_clock::now(); \
      f; \
      t << duration_cast<duration<LongFloat>>(steady_clock::now() - t0).count(); \
    }while(0)
    
    std::cerr << "progress (log2): " << std::flush;
    const uint64_t nSweeps = _n_sweeps + _thermalization;
    int done = ceil(-log(nSweeps)/log(2));
    for ( ; _sweep_num<nSweeps; ++_sweep_num ) {
      _thermalizing = _sweep_num < _thermalization;
      switch (_update_method) {
        case UpdateMethod::local:
          TIME(_t_local, local_update(false));
          if (!_thermalizing)
            TIME(_t_measure, measure());
          break;
        case UpdateMethod::global:
          if (n_ == 1)
            TIME(_t_global, global_update(!_thermalizing));
          else {
            FOR(i,4) {
              if ( bool_dist(random_engine) )
                TIME(_t_local,  local_update());
              else
                TIME(_t_global, global_update());
            }
            if ( !_thermalizing )
              TIME(_t_measure, global_update(true));
          }
          break;
        default:
          Assert_(false, _update_method);
      }
      
      if ( (_sweep_num%16)==0 )
        normalize_spins(); // to prevent roundoff error
      
      if ( log(double(_sweep_num+1)/nSweeps)/log(2) > done ) {
        std::cerr << done << " " << std::flush;
        ++done;
      }
    }
    std::cerr << std::endl;
    #undef TIME
  }
  
  void warn(string warning)
  {
    std::cerr << "WARNING: " << warning << std::endl;
    warnings.push_back('"' + warning + '"');
  }
  
  const Size         N; // number of spins
  const vector<uint> L_;
  const uint         n_;
  Float              J_;
  vector<Float>      J_anisotropy;
  Float              JJ;
  
  friend ostream& operator<<(ostream &os, const MC &mc);
  friend ostream& operator<<(ostream &os, UpdateMethod method);
  
  virtual string showSpins() const =0;
  virtual string layerMeans() const =0;
  
protected:
  MC(Size N_, const vector<uint> &L_v, uint n__)
    : N(N_),
      L_(L_v),
      n_(n__),
      J_(0),
      J_anisotropy(L_v.size(), 1),
      JJ(0),
     _update_method(UpdateMethod::global),
      index_dist(0, N_-1),
     _layer_dims(0),
     _thermalizing(false),
     _thermalization(0),
     _annealing(0),
     _n_sweeps(0),
     _sweep_num(0),
     _sum_1(0,1),
     _sum_1l(0,1),
     _sum_2(0,1),
     _sum_2AF(0,1),
     _sumSS_1(0,1),
     _sumSS_1l(0,1),
     _sumSS_1b(0,1),
     _sumSS_1c(0,1)
  { }
  
  virtual void set_potential_func(SpinFunc *V) =0;
  
  virtual void  local_update(bool allow_simple=true) =0;
  virtual void global_update(bool measure=false) =0;
  virtual void measure() =0;
  
  UpdateMethod _update_method;
  
  std::uniform_int_distribution<Index> index_dist;
  
  uint     _layer_dims;
  string   _potential_name;
  bool     _thermalizing;
  uint64_t _thermalization, _annealing, _n_sweeps, _sweep_num;
  
  Value         _sum_1, _sum_1l, _sum_2, _sum_2AF, _sum_2_q1, _sum_2_q2, _sum_4, _sum_6;
  Value         _sumSS_1, _sumSS_1l, _sumSS_2, _sumSS_1b, _sumSS_2b, _sumSS_1c, _sumSS_2c, _sumSS_2_q1, _sumSS_2_q2, _sumSS_2b_q1l, _sumSS_2b_q2l, _sumSS_2b_q1i, _sumSS_2b_q2i, _sumSS_2c_q1l, _sumSS_2c_q2l, _sumSS_2c_q1i, _sumSS_2c_q2i;
  Value         _n_measurement_clusters, _n_potential_clusters;
  vector<Value> _n_clusters;
  Value         _t_local, _t_global, _t_measure;

  vector<string> warnings;
};
MC::~MC() {}

ostream& operator<<(ostream &os, MC::UpdateMethod method)
{
  String method_str = nullptr;
  switch (method) {
    case MC::UpdateMethod:: local: method_str =  "local"; break;
    case MC::UpdateMethod::global: method_str = "global"; break;
    default:
      Assert_(false, int(method));
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
  
  Assert_(false, method_str);
}

// Spin

template<uint n, typename _float=Float>
struct Spin_
{
  Spin_() = default;
  Spin_(const Spin_ &s) = default;
  Spin_(const array<_float,n> &s) : _s(s) {}
  
  template<typename _float_> explicit
  Spin_(const Spin_<n,_float_> &s)
  { FOR(a,n) _s[a] = s[a]; }
  
  _float& operator[](uint a)       { return _s[a]; }
  _float  operator[](uint a) const { return _s[a]; }
  
  template<typename _float_>
  void  operator+=(Spin_<n,_float_> s)       { FOR(a,n) _s[a] += s[a]; }
  void  operator-=(Spin_            s)       { FOR(a,n) _s[a] -= s[a]; }
  Spin_ operator+ (Spin_            s) const { s += _this; return s; }
  Spin_ operator- (Spin_            s) const { Spin_ ret; FOR(a,n) ret[a] =  _this[a] - s[a]; return ret; }
  Spin_ operator- (                  ) const { Spin_ ret; FOR(a,n) ret[a] = -_this[a]       ; return ret; }
  
  void  operator*=(_float x)       { FOR(a,n) _s[a] *= x; }
  void  operator/=(_float x)       { _this *= 1/x; }
  
  _float operator|(Spin_ s) const { _float x=0; FOR(a,n) x += _s[a]*s[a]; return x; } // dot product
  
  _float norm() const { return sqrt(norm2(_this)); }
  
  void normalize() {
    if (n>1)
      _this /= norm();
    else
      _s[0] = _s[0]>0 ? 1 : -1;
  }
  
  static Spin_ random()
  {
    Spin_ s;
    if ( n > 1 ) {
      Float norm_;
      do {
        FOR(a,n) s[a] = normal_dist(random_engine);
        norm_ = s|s;
      } while (norm_ < .01);
      s /= sqrt(norm_);
    }
    else
      s[0] = 2*bool_dist(random_engine) - 1;
    return s;
  }
  
private:
  template<uint n_, typename _float_>
  friend ostream& operator<<(ostream &os, const Spin_<n_,_float_> &s);
  
  array<_float,n> _s;
};

template<uint n, typename _float>
Spin_<n,_float> operator*(_float x, Spin_<n,_float> s) { s *= x; return s; }

template<uint n, typename _float>
Spin_<n,_float> operator/(Spin_<n,_float> s, _float x) { s *= 1/x; return s; }

template<uint n, typename _float>
_float norm2(Spin_<n,_float> s) { return s|s; }

template<uint n, typename _float>
ostream& operator<<(ostream &os, const Spin_<n,_float> &s)
{ return os << s._s; }

// SpinFlipper

struct SpinFlipper {
  virtual ~SpinFlipper();
  virtual void reset() =0;
  
  virtual uint  n_flippers() const { return 1; }
  virtual uint flipper_num() const { return 0; }
  virtual bool allow_invert() const { return true; } // true if s -> -s is a symmetry
};
SpinFlipper::~SpinFlipper() {}

template<uint n>
struct SpinFlipper_ : public SpinFlipper {
  typedef Spin_<n> Spin;
  virtual ~SpinFlipper_() {}
  
  void flip(Spin &s, bool measure) const { if (measure||n==1) s=-s; else flip_(s); }
  Spin flipped(Spin s, bool measure) const { flip(s,measure); return s; };
  
protected:
  virtual void flip_(Spin &s) const =0; // must satisfy V(R.s) = V(s), and is assumed to be linear
};

template<uint n>
struct InvertSpin_ : public SpinFlipper_<n> {
  typedef Spin_<n> Spin;
  
  virtual void  reset() final { r = Spin::random(); }
  virtual void  flip_(Spin &s) const final { s += (-2*(s|r)) * r; }
  
  Spin r;
};

template<uint n>
struct SpinOperator_ : public SpinFlipper_<n> {
  typedef Spin_<n>              Spin;
  typedef function<void(Spin&)> Operator;
  
  SpinOperator_(const vector<Operator> &&operators_, bool allow_invert_=true)
    : operators(operators_),
     _num(0),
     _dist(0, operators.size()-1),
     _allow_invert(allow_invert_)
  { }
  
  virtual uint  n_flippers() const final { return operators.size(); }
  virtual uint flipper_num() const final { return _num; }
  
  virtual void reset() final { _num = _dist(random_engine); }
  virtual void flip_(Spin &s) const final { operators[_num](s); }
  
  virtual bool allow_invert() const final { return _allow_invert; }
  
  const vector<Operator> operators;
  
private:
  uint _num;
  std::uniform_int_distribution<uint> _dist;
  const bool _allow_invert;
};

// #include "SpinOperatorData.hh"
SpinOperator_<1> ising_flipper({ [](Spin_<1> &s){s = -s;} });

// SpinMatrix

template<uint n>
using SpinMatrix_ = array<Spin_<n>,n>;

template<uint n>
Spin_<n> operator|(const SpinMatrix_<n> &M, const Spin_<n> s0) {
  Spin_<n> s = s0;
  FOR(a,n) s[a] = M[a] | s0;
  return s;
}

template<uint n>
SpinMatrix_<n> operator|(const SpinMatrix_<n> &M1, const SpinMatrix_<n> &M2)
{
  SpinMatrix_<n> M;
  
  FOR(a,n)
  FOR(b,n) {
    LongFloat m = 0;
    FOR(c,n)
      m += M1[a][c] * M2[c][b];
    M[a][b] = m;
  }
  
  return M;
}

template<uint n>
LongFloat tr(const SpinMatrix_<n> &M) {
  LongFloat x = 0;
  FOR(a,n)
    x += M[a][a];
  return x;
}

template<uint n>
LongFloat tr(const SpinMatrix_<n> &M1, const SpinMatrix_<n> &M2) {
  LongFloat x = 0;
  FOR(a,n)
  FOR(b,n)
    x += M1[a][b] * M2[b][a];
  return x;
}

// SpinFunc

struct SpinFunc {
  virtual ~SpinFunc();
  virtual vector<Float> coefficients() const =0;
  
  bool nonzero() const {
    vector<Float> cs = coefficients();
    for(Float c: cs)
      if ( c != 0 )
        return true;
    return false;
  }
  
  vector<string> measurement_names;
  vector<Value>  measurements;
};
SpinFunc::~SpinFunc() {}

const bool measureSSQ = true;

template<uint n>
struct SpinFunc_ : public SpinFunc {
  typedef Spin_<n,    Float> Spin;
  typedef Spin_<n,LongFloat> SpinSum;
  
  virtual ~SpinFunc_() {}
  virtual Float operator()(Spin s) const =0;
  virtual Spin ideal_spin() const = 0;
  
  struct MeasureArgs {SpinSum sum; LongFloat sum_4; const Spin *spins; MC::Size N;};
  virtual void measure(MeasureArgs r) {}
  
  Spin ideal_spin_from(vector<Spin> &candidates) const
  {
    Spin *min_s = nullptr;
    Float min_V = infinity<Float>();
    for (Spin &s0: candidates) {
      s0.normalize();
      const Float V0 = _this(s0);
      if ( V0 < min_V ) {
        min_s = &s0;
        min_V = V0;
      }
    }
    Assert (min_s);
    return *min_s;
  }
};

// MC_

template<uint dim, // # of spacial dimensions
         uint n>   // # of spin components
class MC_ : public MC
{
  enum class InitTemp { none, high, low };
public:
  typedef Spin_<n,    Float>    Spin;
  typedef Spin_<n,LongFloat>    SpinSum;
  
  // POD is zero initialized by {}
  static_assert( std::is_pod<Spin>::value, "Spin isn't POD" );
  
  // lattice
  
  typedef array<uint,dim> Pos;
  
  struct Neighbor { Index i; Pos p; uint d; };
  typedef vector<Neighbor> NearestNeighbors;
  
  Index index(const Pos p) const // todo: consider blocking
  {
    Index i = p[0];
    for (uint d=1; d<dim; ++d)
      i = i*L[d] + p[d];
    return i;
  }
  
  Pos pos(Index i) const // todo: consider blocking
  {
    Pos p;
    for (uint d=dim-1; d>0; --d) {
      p[d] = i%L[d];
      i /= L[d];
    }
    p[0] = i;
    return p;
  }
  
  bool AF_signQ(Pos p) const
  {
    uint tot = 0;
    for (uint d=_layer_dims; d<dim; ++d)
      tot += p[d];
    return tot % 2;
  }
  
  bool AF_layer_signQ(Pos p) const
  {
    uint tot = 0;
    FOR(d, _layer_dims)
      tot += p[d];
    return tot % 2;
  }
  
  uint mod_d(uint x, uint d) const { return (x + L[d]) % L[d]; }
  
  NearestNeighbors nearestNeighbors(const Pos p) const
  {
    NearestNeighbors nn;
    
    for (uint d=_layer_dims; d<dim; ++d)
    for (int dir=-1; dir<=+1; dir+=2) {
      Pos q = p;
      q[d] = mod_d(q[d] + dir, d);
      const Index j = index(q);
      
      nn.push_back(Neighbor{j, q, d});
    }
    
    return nn;
  }
  
  // MC
  
  virtual void clear_spins() final {
    _init_temp = InitTemp::low;
    Spin s;
    if ( _V )
      s = _V->ideal_spin();
    else
      FOR(a,n) s[a] = (a==0);
    FOR(i,N) _spins[i] = s;
  }
  
  virtual void randomize_spins() final
  { _init_temp = InitTemp::high; FOR(i,N) _spins[i] = Spin::random(); }
  
  virtual void normalize_spins() final
  { FOR(i,N) _spins[i].normalize(); }
  
  Float Jnn(Pos ip, Neighbor nn) const
  {
    Float J0 = 0;
    for (uint d=0; d<_layer_dims; ++d)
    for (int dir=-1; dir<=+1; dir+=2) {
      Pos p1=ip, p2=nn.p;
      p1[d] = mod_d(p1[d] + dir, d);
      p2[d] = mod_d(p2[d] + dir, d);
      J0 += _spins[index(p1)] | _spins[index(p2)];
    }
    return J[nn.d] + JJ*J0;
  }
  
  virtual void local_update(const bool allow_simple) final __attribute__((hot))
  {
    if ( _setup )
      setup();
    
    const Size n_updates = 1+index_dist(random_engine);
    for (Size count=0; count<n_updates; ++count) {
      const Index i      = index_dist(random_engine);
      const bool  simple = n>1 && allow_simple && bool_dist(random_engine);
      
      const Spin s1 = _spins[i];
      Spin sum_nn{};
      const Pos ip = pos(i);
      for (Neighbor nn : nearestNeighbors(ip))
        sum_nn += Jnn(ip,nn) * _spins[nn.i];
      
      Spin s2;
      Float delta_E = 0;
      if (!simple) {
        s2 = n==1 ? -s1 : Spin::random();
        delta_E += -((s2-s1)|sum_nn);
      } else {
        const Float norm = norm2(sum_nn);
        s2 = norm!=0 ? (2*(s1|sum_nn)/norm)*sum_nn - s1 : -s1;
      }
      
      if (_V)
        delta_E += V(s2) - V(s1);
      
      delta_E = anneal(delta_E);
      if ( delta_E <= 0 || uniform_dist(random_engine) < exp(-delta_E) )
        _spins[i] = s2;
    }
  }
  
  virtual void global_update(const bool measure_=false) final __attribute__((hot))
  {
    if ( _setup )
      setup();
    
    if ( _setup_global_update )
    {
      FOR(d, dim)
      if ( _use_q_dim[d] )
        _clusterSums_q[d].reset(new array<SpinSum,4>[N]); 
      
      _n_clusters.resize( _flipper->n_flippers() );
      check_potential();
      _setup_global_update = false;
    }
    
    if ( measure_ && !_flipper->allow_invert() ) {
      // s -> -s improved estimator code isn't allowed
      Assert(n>1);
      measure();
      return;
    }
    
    // use _wolff_flipper despite the symmetry breaking potential
    const bool force_wolff = !measure_ && _V && bool_dist(random_engine);
    SpinFlipper_<n> &flipper = force_wolff ? _wolff_flipper : *_flipper;
    
    if ( !measure_ )
      flipper.reset();
    _cluster.assign(N, false);
    
    Size nClusters = 0;
    FOR(i0, N)
    if ( !_cluster[i0] ) {
      const Pos ip0 = pos(i0);
      _clusterLayers[nClusters] = ip0;
      const bool flip_cluster = !force_wolff && bool_dist(random_engine);
      Index *newIndex = _newIndexStack.get();
      
      #define sum_cluster_q_(op, i, p) \
        do { \
          _clusterSums  [nClusters] op (s); \
          _clusterSumsAF[nClusters] op (AF_signQ(p) ? s : -s); \
           \
          FOR(d, dim) \
          if ( _use_q_dim[d] ) { \
            const uint x1 =    p[d]; \
            const uint x2 = (2*p[d])%L[d]; \
            _clusterSums_q[d][nClusters][0] op (_cos[x1] * s); \
            _clusterSums_q[d][nClusters][1] op (_sin[x1] * s); \
            _clusterSums_q[d][nClusters][2] op (_cos[x2] * s); \
            _clusterSums_q[d][nClusters][3] op (_sin[x2] * s); \
          } \
        } while(0)
      
      #define set_cluster_q(i, p) sum_cluster_q_(=SpinSum, i, p);
      #define add_cluster_q(i, p) sum_cluster_q_(+=, i, p);
      
      *newIndex = i0;
      _cluster[i0] = true;
      uint n_cluster_indices = 0;
      Float  cluster_delta_V = 0;
      
      {
        const Spin s = _spins[i0];
        if ( measure_ )
          set_cluster_q(i0, ip0);
        if ( force_wolff ) {
          _cluster_indices[n_cluster_indices++] = i0; // _cluster_indices is a list of indices in the cluster
          cluster_delta_V += V(_wolff_flipper.flipped(s,false)) - V(s);
        }
        if ( flip_cluster )
          flipper.flip(_spins[i0], measure_);
      }
      
      do {
        const Index j  = *(newIndex--);
        const Pos   jp = pos(j);
        
        for(Neighbor nn : nearestNeighbors(jp)) {
          const Index i  = nn.i;
          const Pos   ip = pos(i);
          if ( !_cluster[i] ) {
            const Spin s = _spins[i];
            const Spin flipped_s = flipper.flipped(s, measure_);
            const Float delta_E = anneal( -Jnn(jp,nn)*(1-2*flip_cluster) * ((flipped_s-s)|_spins[j]) );
            if ( delta_E > 0 && uniform_dist(random_engine) > exp(-delta_E) ) {
              *(++newIndex) = i;
              _cluster[i] = true;
              if ( measure_ )
                add_cluster_q(i, ip);
              if ( force_wolff ) {
                _cluster_indices[n_cluster_indices++] = i;
                cluster_delta_V += V(flipped_s) - V(s);
              }
              if ( flip_cluster )
                _spins[i] = flipped_s;
            }
          }
        }
      }
      while ( newIndex+1 != _newIndexStack.get() );
      
      #undef sum_cluster_q_
      #undef set_cluster_q
      #undef add_cluster_q
      
      if ( force_wolff ) {
        cluster_delta_V = anneal(cluster_delta_V);
        if ( cluster_delta_V <= 0 || uniform_dist(random_engine) < exp(-cluster_delta_V) )
          FOR(k, n_cluster_indices)
            _wolff_flipper.flip(_spins[_cluster_indices[k]],false);
      }
      
      ++nClusters;
    }
    
    if ( measure_ )
    {
      check_thermalization();
      
      typedef SpinMatrix_<n> SpinMatrix;
      static_assert( std::is_pod<SpinMatrix>::value, "SpinMatrix isn't POD" );
      
      LongFloat sum_2l = 0;
      for (uint c=0; c<nClusters; )
      {
        LongFloat  S6=0, S4_1=0, S2AF=0;
        SpinSum    S1{}, S1AF{}, Q2{};
        SpinMatrix M2{}, M4{};
        LongFloat  sum_2_q1=0, sum_2_q2=0;
        const Pos p0 = _clusterLayers[c];
        while (true) {
          const SpinSum   s   = _clusterSums[c],
                          sAF = _clusterSumsAF[c];
          const LongFloat s2  = s|s;
          
          S1   += s;
          S1AF += sAF;
          S2AF += sAF|sAF;
          S6 += Pow<3>(s2);
          FOR(a, n) {
            if ( n > 1 ) {
              S4_1  += Pow<4>(s[a]);
              Q2[a] += Pow<2>(s[a]);
            }
            FOR(b, n) {
              M2[a][b] +=      s[a]*s[b];
              M4[a][b] += s2 * s[a]*s[b];
            }
          }
          
          FOR(d, dim)
          if ( _use_q_dim[d] )
          FOR(cs, 2) {
            sum_2_q1 += norm2(_clusterSums_q[d][c][0+cs]);
            sum_2_q2 += norm2(_clusterSums_q[d][c][2+cs]);
          }
          
          ++c;
          bool done = c>=nClusters;
          FOR(d, _layer_dims)
            done = done || p0[d] != _clusterLayers[c][d];
          if (done)
            break;
        }
        const SpinMatrix M2_2 = M2|M2;
        
        LongFloat N1 = 1;
        for (uint d=_layer_dims; d<dim; ++d)
          N1 *= L[d];
        S1   /= N1;
        S1AF /= N1;
        const LongFloat
        N2      = sq(N1),
        N4      = sq(N2),
        N6      = N2*N4,
        sum_2   = norm2(S1),
        sum_2AF = norm2(S1AF),
        sum_1   = sqrt(sum_2),
        S2      = tr(M2),
        S4      = tr(M4),
        S_2_2   = tr(M2_2),
        S_2_2_2 = tr(M2_2, M2),
        S_4_2   = tr(M4, M2);
        
        LongFloat sum_4;
        _sum_1   .add       (sum_1);
        sum_2l += sum_2;
        _sum_2   .add       (S2  /N2, sum_2);
        _sum_2AF .add       (S2AF/N2, sum_2AF);
        _sum_2_q1 <<         sum_2_q1 / (_n_q_dim*N2);
        _sum_2_q2 <<         sum_2_q2 / (_n_q_dim*N2);
        _sum_4    << (sum_4=(sq(S2) - 2*S4 + 2*S_2_2)/N4 );
        _sum_6    <<        (Pow<3>(S2) - 6*S2*S4 + 16*S6 + 6*S2*S_2_2 - 24*S_4_2 + 8*S_2_2_2)/N6;
        if (_V)
          _V->measure({S1, sum_4, _spins.get(), N});
      }
      
      uint Nl = 1;
      FOR(d, _layer_dims)
        Nl *= L[d];
      _sum_1l.add(sqrt(sum_2l/Nl));
      
      _n_measurement_clusters << nClusters;
      
      measureSS();
    }
    else if ( !_thermalizing ) {
      if ( force_wolff )
        _n_potential_clusters << nClusters;
      else
        _n_clusters[_flipper->flipper_num()] << nClusters;
    }
  }
  
  virtual void set_flipper(SpinFlipper *flipper) final
  {
    if (flipper) {
      _flipper = dynamic_cast<SpinFlipper_<n>*>(flipper);
      Assert_(_flipper, n);
    }
    else
      _flipper = n>1 ? &_wolff_flipper : dynamic_cast<SpinFlipper_<n>*>(n==1 ? &ising_flipper : nullptr);
  }
  
  virtual const SpinFunc* potential() const final { return _V; }
  
  virtual string showSpins() const final
  {
    stringstream ret;
    ret << "{";
    FOR(i, N) {
      if (n==1)
        ret << _spins[i][0];
      else
        ret << _spins[i];
      if ( i+1 < N )
        ret << ", ";
    }
    ret << "}\n";
    return ret.str();
  }
  
  virtual string layerMeans() const final
  {
    uint Ni = 1;
    for (uint d=_layer_dims; d<dim; ++d)
      Ni *= L[d];
    
    stringstream ret;
    ret << "{";
    SpinSum sum{};
    FOR(i, N) {
      sum += _spins[i];
      
      if ( (i+1)%Ni == 0 ) {
        sum /= Ni;
        if (n==1)
          ret << sum[0];
        else
          ret << sum;
        sum = SpinSum();
        
        if ( i+1 < N )
          ret << ", ";
      }
    }
    ret << "}\n";
    return ret.str();
  }
  
  Float V(Spin s) const { return (*_V)(s); }
  
protected:
  void measureSS()
  {
    if (!(0 < _layer_dims &&_layer_dims < dim && measureSSQ))
      return;
    
    uint Ni = 1;
    for (uint d=_layer_dims; d<dim; ++d)
      Ni *= L[d];
    uint    n_sumSS_2l=0;
    LongFloat sumSS_2l=0, sumSSb[dim]{}, sumSS[dim]{}, sumSSb_q1[dim][dim][2]{}, sumSSb_q2[dim][dim][2]{}, sumSS_q1[dim][dim][2]{}, sumSS_q2[dim][dim][2]{};
    FOR(i, N) {
      const Pos p = pos(i);
      const Spin s = _spins[i];
      const Float sgn = AF_layer_signQ(p) ? 1 : -1;
      
      FOR(d0, dim)
      if ((d0 < _layer_dims && L[d0] == L[0]) || _use_q_dim[d0]) {
        Pos p2 = p;
        p2[d0] = mod_d(p2[d0] + 1, d0);
        const Spin s2 = _spins[index(p2)];
        
        const Float ss  = s|s2;
        const Float ssb = sgn * ss;
        sumSSb[d0] += ssb;
        sumSS [d0] += ss;
        
        FOR(d, dim)
        if ( (d0 < _layer_dims && _use_q_dim[d]) || (d0 >= _layer_dims && d != d0) ) {
          const uint x1 =    p[d];
          const uint x2 = (2*p[d])%L[d];
          const Float *cos_ = d < _layer_dims ? _cos0.get() : _cos.get(),
                      *sin_ = d < _layer_dims ? _sin0.get() : _sin.get();
          sumSSb_q1[d0][d][0] += cos_[x1] * ssb;
          sumSSb_q1[d0][d][1] += sin_[x1] * ssb;
          sumSSb_q2[d0][d][0] += cos_[x2] * ssb;
          sumSSb_q2[d0][d][1] += sin_[x2] * ssb;
          sumSS_q1 [d0][d][0] += cos_[x1] * ss;
          sumSS_q1 [d0][d][1] += sin_[x1] * ss;
          sumSS_q2 [d0][d][0] += cos_[x2] * ss;
          sumSS_q2 [d0][d][1] += sin_[x2] * ss;
        }
      }
      
      if ( (i+1)%Ni == 0 )
      FOR(d0, _layer_dims)
      if (L[d0] == L[0]) {
        LongFloat sumSS_2_q1=0, sumSS_2_q2=0;
        uint nSS_q_dim=0;
        for (uint d=_layer_dims; d<dim; ++d)
        if (_use_q_dim[d]) {
          ++nSS_q_dim;
          FOR(cs, 2) {
            sumSS_2_q1 += sq(sumSS_q1[d0][d][cs]);
                             sumSS_q1[d0][d][cs] = 0;
            sumSS_2_q2 += sq(sumSS_q2[d0][d][cs]);
                             sumSS_q2[d0][d][cs] = 0;
          }
        }
        const LongFloat Ni2 = sq(LongFloat(Ni)),
                    sumSS_2 = sq(sumSS[d0]/Ni);
                                 sumSS[d0] = 0;
        
        _sumSS_1   .add(sqrt(sumSS_2));
        ++n_sumSS_2l;
        sumSS_2l += sumSS_2;
        _sumSS_2    <<  sumSS_2;
        _sumSS_2_q1 <<  sumSS_2_q1 / (nSS_q_dim*Ni2);
        _sumSS_2_q2 <<  sumSS_2_q2 / (nSS_q_dim*Ni2);
      }
    }
    
    _sumSS_1l.add(sqrt(sumSS_2l/n_sumSS_2l));
    
    for (uint d0=_layer_dims; d0<dim; ++d0)
    if (_use_q_dim[d0]) {
      LongFloat sumSS_2b_q1l=0, sumSS_2b_q2l=0, sumSS_2c_q1l=0, sumSS_2c_q2l=0;
      uint nSSb_ql_dim=0;
      FOR(d, _layer_dims)
      if (L[d] == L[0]) {
        ++nSSb_ql_dim;
        FOR(cs, 2) {
          sumSS_2b_q1l += sq(sumSSb_q1[d0][d][cs]);
          sumSS_2b_q2l += sq(sumSSb_q2[d0][d][cs]);
          sumSS_2c_q1l += sq( sumSS_q1[d0][d][cs]);
          sumSS_2c_q2l += sq( sumSS_q2[d0][d][cs]);
        }
      }
      
      LongFloat sumSS_2b_q1i=0, sumSS_2b_q2i=0, sumSS_2c_q1i=0, sumSS_2c_q2i=0;
      uint nSSb_qi_dim=0;
      for (uint d=_layer_dims; d<dim; ++d)
      if (d != d0 && L[d] == L[_layer_dims]) {
        ++nSSb_qi_dim;
        FOR(cs, 2) {
          sumSS_2b_q1i += sq(sumSSb_q1[d0][d][cs]);
          sumSS_2b_q2i += sq(sumSSb_q2[d0][d][cs]);
          sumSS_2c_q1i += sq( sumSS_q1[d0][d][cs]);
          sumSS_2c_q2i += sq( sumSS_q2[d0][d][cs]);
        }
      }
      
      const LongFloat N2 = sq(LongFloat(N)),
                sumSS_2b = sq(sumSSb[d0]/N),
                sumSS_2c = sq( sumSS[d0]/N);
      
      _sumSS_1b.add(sqrt(sumSS_2b));
      _sumSS_2b      <<  sumSS_2b;
      _sumSS_2b_q1l  <<  sumSS_2b_q1l / (nSSb_ql_dim*N2);
      _sumSS_2b_q2l  <<  sumSS_2b_q2l / (nSSb_ql_dim*N2);
      _sumSS_2b_q1i  <<  sumSS_2b_q1i / (nSSb_qi_dim*N2);
      _sumSS_2b_q2i  <<  sumSS_2b_q2i / (nSSb_qi_dim*N2);
      
      _sumSS_1c.add(sqrt(sumSS_2c));
      _sumSS_2c      <<  sumSS_2c;
      _sumSS_2c_q1l  <<  sumSS_2c_q1l / (nSSb_ql_dim*N2);
      _sumSS_2c_q2l  <<  sumSS_2c_q2l / (nSSb_ql_dim*N2);
      _sumSS_2c_q1i  <<  sumSS_2c_q1i / (nSSb_qi_dim*N2);
      _sumSS_2c_q2i  <<  sumSS_2c_q2i / (nSSb_qi_dim*N2);
    }
  }
  
  virtual void measure() final
  {
    check_thermalization();
    
    SpinSum sum{}, sumAF{}, sum_q1[dim][2]{}, sum_q2[dim][2]{};
    FOR(i, N) {
      const Pos p = pos(i);
      bool sumQ = true;
      FOR(d, _layer_dims)
        sumQ = sumQ && p[d]==0;
      if ( sumQ ) { // this could be factored into the FOR(i, N) loop
        const Spin s = _spins[i];
        sum += s;
        sumAF += AF_signQ(p) ? s : -s;
        
        FOR(d, dim)
        if ( _use_q_dim[d] ) {
          const uint x1 =    p[d];
          const uint x2 = (2*p[d])%L[d];
          sum_q1[d][0] += _cos[x1] * s;
          sum_q1[d][1] += _sin[x1] * s;
          sum_q2[d][0] += _cos[x2] * s;
          sum_q2[d][1] += _sin[x2] * s;
        }
      }
    }
    LongFloat N1 = 1;
    for (uint d=_layer_dims; d<dim; ++d)
      N1 *= L[d];
    sum /= N1;
    sumAF /= N1;
    
    const LongFloat
    N2 = sq(N1),
    sum_2   = norm2(sum),
    sum_2AF = norm2(sumAF),
    sum_4   =    sq(sum_2);
    
    LongFloat sum_2_q1=0, sum_2_q2=0;
    FOR(d, dim)
      if ( _use_q_dim[d] )
        FOR(cs, 2) {
          sum_2_q1 += norm2(sum_q1[d][cs]);
          sum_2_q2 += norm2(sum_q2[d][cs]);
        }
    sum_2_q1 /= _n_q_dim*N2;
    sum_2_q2 /= _n_q_dim*N2;
    
    const LongFloat sum_1 = sqrt(sum_2);
    
    _sum_1   .add       (sum_1);
    _sum_2   .add       (sum_2);
    _sum_2AF .add       (sum_2AF);
    _sum_2_q1 <<         sum_2_q1;
    _sum_2_q2 <<         sum_2_q2;
    _sum_4    <<         sum_4;
    _sum_6    <<  Pow<3>(sum_2);
    if (_V)
      _V->measure({sum, sum_4, _spins.get(), N});
    
    measureSS();
  }
  
  Float anneal(Float delta_E) const {
    if ( _sweep_num >= _annealing )
      return delta_E;
    else if ( _init_temp == InitTemp::high )
      return (delta_E*_sweep_num)/_annealing;
    else if ( _init_temp == InitTemp::low )
      return (delta_E*_annealing)/(_sweep_num+1);
    else
      Assert(false);
  }
  
  void setup()
  {
    _n_q_dim = 0;
    FOR(d, dim) {
      if ( d < _layer_dims )
        Assert(J_anisotropy[d] == 0);
      J[d] = J_ * J_anisotropy[d];
      _use_q_dim[d] = d >= _layer_dims && L[d] == L[_layer_dims] && J[d] == J[_layer_dims];
      if ( _use_q_dim[d] )
        ++_n_q_dim;
    }
    
    _cos.reset(new Float[L[_layer_dims]]);
    _sin.reset(new Float[L[_layer_dims]]);
    FOR(x, L[_layer_dims]) {
      _cos[x] = cos((2*pi*x)/L[_layer_dims]);
      _sin[x] = sin((2*pi*x)/L[_layer_dims]);
    }
    
    _cos0.reset(new Float[L[0]]);
    _sin0.reset(new Float[L[0]]);
    FOR(x, L[0]) {
      _cos0[x] = cos((2*pi*x)/L[0]);
      _sin0[x] = sin((2*pi*x)/L[0]);
    }
    
    _setup = false;
  }
  
  void check_potential() const
  {
    if ( _V ) {
      const Spin ideal_spin = _V->ideal_spin();
      Assert_( is_small(abs(ideal_spin.norm()-1)), ideal_spin);
      const Float V_min = V(ideal_spin);
      FOR(i, 100) {
        const bool measure0 = i<3 && _flipper->allow_invert();
        _flipper->reset();
        const Spin  s1 = Spin::random();
        const Spin  s2 = _flipper->flipped(s1,measure0);
        const Spin  s3 = _flipper->flipped(s2,measure0); // should == s1
        const Float V1 = V(s1);
        const Float delta_V = V(s2) - V1;
        const Float norm_s13 = (s1-s3).norm();
        Assert_( is_small(delta_V), s1, s2, delta_V );
        Assert_( is_small(norm_s13), s1, s3, norm_s13 );
        Assert_( V_min <= V1, V_min, V1, ideal_spin, s1 );
      }
    }
  }
  
  void check_thermalization()
  {
    if ( _check_thermalization ) {
      _check_thermalization = false;
      uint count = 0;
      if ( _init_temp == InitTemp::low && _V && _V->nonzero() ) {
        const Float V_min = V(_V->ideal_spin()) + sqrt(std::numeric_limits<Float>::epsilon());
        FOR(i, N)
          count += V(_spins[i]) <= V_min;
      }
      if ( count )
        warn("un-thermalized spins: " + stringify(count));
    }
  }
  
  MC_(const array<uint,dim> &L__, const vector<uint> &L_v, Size N_)
    : MC(N_, L_v, n),
      L(L__),
     _spins(new Spin[N_]),
     _newIndexStack(new Index[N_]),
     _clusterLayers(new Pos[N_]),
     _clusterSums  (new SpinSum[N_]),
     _clusterSumsAF(new SpinSum[N_]),
     _n_q_dim(0),
     _use_q_dim{},
     _cluster_indices(new Index[N_]),
     _flipper(nullptr),
     _V(nullptr),
     _init_temp(InitTemp::none),
     _setup(true),
     _setup_global_update(true),
     _check_thermalization(true)
  {
    set_flipper(nullptr);
    
    FOR(d, dim)
      J[d] = NAN;
  }
  
  virtual void set_potential_func(SpinFunc *new_V) final
  {
    _V = dynamic_cast<SpinFunc_<n>*>(new_V);
    
    if (new_V)
      Assert_(_V, n);
  }
  
  friend class MC;
  
  virtual ~MC_() {}
  
  const array<uint,dim> L; // lengths
  array<Float,dim>      J; // includes anisotropy
  
private:
  const unique_ptr<Spin[]>       _spins;
  vector<bool>                   _cluster;
  const unique_ptr<Index[]>      _newIndexStack;
  const unique_ptr<Pos[]>        _clusterLayers;
  const unique_ptr<SpinSum[]>    _clusterSums, _clusterSumsAF;
  unique_ptr<array<SpinSum,4>[]> _clusterSums_q[dim]; // 4 = #(q1, q2) * #(cos, sin)
  uint                           _n_q_dim;
  bool                           _use_q_dim[dim];
  unique_ptr<Index[]>            _cluster_indices;
  unique_ptr<Float[]>            _cos, _sin, _cos0, _sin0;
  
  InvertSpin_<n>      _wolff_flipper;
  SpinFlipper_<n>    *_flipper;
  SpinFunc_<n>       *_V;
  
  InitTemp _init_temp;
  bool _setup, _setup_global_update, _check_thermalization;
};

bool __brief_values = false;
bool __print_spins  = false;
ostream& operator<<(ostream &os, const MC &mc)
{
  const SpinFunc *const potential = mc.potential();
  vector<IOFloat> J{mc.J_}, J_anisotropy;
  if (potential) {
    const vector<Float> coefficients = potential->coefficients();
    J.insert(J.end(), coefficients.begin(), coefficients.end());
  }
  for (Float J_a: mc.J_anisotropy)
    J_anisotropy.push_back(J_a);
  IOFloat JJ = mc.JJ;
  
  const int old_prec = os.precision();
  os.precision(15);
  Value::brief = __brief_values;
  os << "{\n"
        "\"L\" -> " << mc.L_ << ",\n"
        "\"n\" -> " << mc.n_ << ",\n"
        "\"J\" -> " <<    J  << ",\n"
        "\"J anisotropy"      "\" -> "   <<   J_anisotropy              <<   ",\n"
        "\"JJ"                "\" -> "   <<   JJ                        <<   ",\n"
        "\"layer dims"        "\" -> "   << mc._layer_dims              <<   ",\n"
        "\"potential"         "\" -> \"" << mc._potential_name          << "\",\n"
        "\"update method"     "\" -> \"" << mc._update_method           << "\",\n"
        "\"sweeps"            "\" -> "   << mc._n_sweeps                <<   ",\n"
        "\"thermalization"    "\" -> "   << mc._thermalization          <<   ",\n"
        "\"annealing"         "\" -> "   << mc._annealing               <<   ",\n"
        "\"|S|"               "\" -> "   << mc._sum_1                   <<   ",\n"
        "\"|S|l"              "\" -> "   << mc._sum_1l                  <<   ",\n"
        "\"S^2"               "\" -> "   << mc._sum_2                   <<   ",\n";
  if ( mc._layer_dims < mc.L_.size() )
  os << "\"S^2_q1"            "\" -> "   << mc._sum_2_q1                <<   ",\n"
        "\"S^2_q2"            "\" -> "   << mc._sum_2_q2                <<   ",\n";
  os << "\"S^2 AF"            "\" -> "   << mc._sum_2AF                 <<   ",\n";
  if ( 0 < mc._layer_dims && mc._layer_dims < mc.L_.size() && measureSSQ ) {
  os << "\"|SS|"              "\" -> "   << mc._sumSS_1                 <<   ",\n"
        "\"|SS|l"             "\" -> "   << mc._sumSS_1l                <<   ",\n"  // rms over layers
        "\"SS^2"              "\" -> "   << mc._sumSS_2                 <<   ",\n"  // (sum_i s_li * s_l'i)^2
        "\"SS^2_q1"           "\" -> "   << mc._sumSS_2_q1              <<   ",\n"  // vs i
        "\"SS^2_q2"           "\" -> "   << mc._sumSS_2_q2              <<   ",\n"
        "\"|SSb|"             "\" -> "   << mc._sumSS_1b                <<   ",\n"
        "\"SS^2b"             "\" -> "   << mc._sumSS_2b                <<   ",\n"  // (sum_li (-)^l s_li * s_li')^2 // intermediate order param
        "\"SS^2b_q1l"         "\" -> "   << mc._sumSS_2b_q1l            <<   ",\n"  // vs l
        "\"SS^2b_q2l"         "\" -> "   << mc._sumSS_2b_q2l            <<   ",\n"
        "\"|SSc|"             "\" -> "   << mc._sumSS_1c                <<   ",\n"
        "\"SS^2c"             "\" -> "   << mc._sumSS_2c                <<   ",\n"  // (s_li * s_li')^2 // intermediate order param
        "\"SS^2c_q1l"         "\" -> "   << mc._sumSS_2c_q1l            <<   ",\n"  // vs l
        "\"SS^2c_q2l"         "\" -> "   << mc._sumSS_2c_q2l            <<   ",\n";
  if ( mc.L_.size() - mc._layer_dims >= 2 )
  os << "\"SS^2b_q1i"         "\" -> "   << mc._sumSS_2b_q1i            <<   ",\n"  // vs i''
        "\"SS^2b_q2i"         "\" -> "   << mc._sumSS_2b_q2i            <<   ",\n"
        "\"SS^2c_q1i"         "\" -> "   << mc._sumSS_2c_q1i            <<   ",\n"  // vs i''
        "\"SS^2c_q2i"         "\" -> "   << mc._sumSS_2c_q2i            <<   ",\n"; }
  os << "\"S^4"               "\" -> "   << mc._sum_4                   <<   ",\n"
        "\"S^6"               "\" -> "   << mc._sum_6                   <<   ",\n";
  
  if (__print_spins)
    os << "\"spins\" -> " << mc.showSpins() << ",\n";
  
  if (0 < mc._layer_dims)
    os << "\"layer means\" -> " << mc.layerMeans() << ",\n";
  
  if (potential) {
    const uint n = potential->measurements.size();
    Assert( potential->measurement_names.size() == n );
    FOR(i, n)
      os << '"' << potential->measurement_names[i] << "\" -> " << potential->measurements[i] << ",\n";
  }
  
  if ( mc.warnings.size() )
    os << "\"warnings\" -> " << mc.warnings << ",\n";
  
  Value::brief = true;
  os << "\"measurement clusters\" -> "   << mc._n_measurement_clusters  <<   ",\n"
        "\"potential clusters""\" -> "   << mc._n_potential_clusters    <<   ",\n"
        "\"clusters"          "\" -> "   << mc._n_clusters              <<   ",\n"
        "\"t local"           "\" -> "   << mc._t_local                 <<   ",\n"
        "\"t global"          "\" -> "   << mc._t_global                <<   ",\n"
        "\"t measure"         "\" -> "   << mc._t_measure               <<    "\n";
  Value::brief = false;
  os.precision(old_prec);
  return os << "}\n";
}

unique_ptr<MC> MC::make(const vector<uint> &L, uint n_fields)
{
  MC::Size N = 1;
  const uint dim = L.size();
  FOR(d,dim) {
    const MC::Size old_N = N;
    N *= L[d];
    Assert_( L[d] && N/L[d] == old_N, old_N, L[d], N ); // check overflow
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
  ELSE_TRY_MC(1,1)
//   ELSE_TRY_MC(1,2)
//ELSE_TRY_MC(1,3)
//ELSE_TRY_MC(1,6)
  ELSE_TRY_MC(2,1)
//ELSE_TRY_MC(2,2)
//ELSE_TRY_MC(2,3)
//ELSE_TRY_MC(2,4)
//ELSE_TRY_MC(2,6)
  ELSE_TRY_MC(3,1)
//ELSE_TRY_MC(3,2)
//ELSE_TRY_MC(3,3)
//ELSE_TRY_MC(3,4)
//ELSE_TRY_MC(3,6)
  ELSE_TRY_MC(4,1)
//ELSE_TRY_MC(4,2)
//ELSE_TRY_MC(4,3)
  else
    Assert_(false, L, N);
  
  #undef ELSE_TRY_MC
  
  return unique_ptr<MC>(mc);
}

// main

int main(const int argc, char *argv[])
{
  using std::cout;
  using std::cerr;
  using std::endl;
  
  // read program options
  namespace po = boost::program_options;
  
  po::options_description generic_options("generic options");
  generic_options.add_options()
      ("help,h",       "print help message")
      ("verbose",      "print verbose output")
      ("brief-values", "don't output detailed measurement errors")
      ("print-spins",  "output final spin configuration");
  
  vector<uint> L;
  uint n;
  vector<IOFloat> J_IO, J_anisotropy_IO;
  IOFloat JJ_IO;
  string potential_name;
  po::options_description system_options("physics options");
  system_options.add_options()
      ("L",            po::value<vector<uint>>(&L)->multitoken(),                  "lengths")
      ("n",            po::value<uint>(&n)->default_value(1),                      "for an O(n) model")
      ("J",            po::value<vector<IOFloat>>(&J_IO)->multitoken(),            "spin coupling and other potential coefficients")
      ("J-anisotropy", po::value<vector<IOFloat>>(&J_anisotropy_IO)->multitoken(), "spin coupling anisotropy for each dimension")
      ("JJ",           po::value<IOFloat>(&JJ_IO)->default_value(0),               "4-spin layer-layer coupling");
//    ("potential",    po::value<string>(&potential_name),                         "potential term: s^4, cos#, OnZ2Z2, vison hexagon|square|triangle c|s-VBS");
  
  uint     n_layer_dims;
  uint64_t n_sweeps=0;
  uint64_t thermalization=0, annealing=0;
  string update_method, initial_state, file_name, replace;
  po::options_description simulation_options("simulation options");
  simulation_options.add_options()
      ("sweeps",         po::value<uint64_t>(&n_sweeps),                               "# of MC sweeps (rounded to a power of 2)")
      ("thermalization", po::value<uint64_t>(&thermalization),                         "# of thermalization sweeps (defaults to #sweeps/4)")
      ("annealing",      po::value<uint64_t>(&annealing),                              "# of annealing sweeps (defaults to #thermalization/4)")
      ("update-method",  po::value<string>(&update_method)->default_value("global"),   "update type: local or global")
      ("initial-state",  po::value<string>(&initial_state)->default_value("low-temp"), "low-temp or high-temp")
      ("layer-dims",     po::value<uint>(&n_layer_dims)->default_value(0),             "number of layer dimensions")
      ("file",           po::value<string>(&file_name),                                "save file (or into directory/)")
      ("replace",        po::value<string>(&replace)->default_value("no"),             "replace file? yes, no, auto, or auto-error");
  
  po::options_description cmdline_options;
  cmdline_options.add(generic_options)
                 .add(system_options)
                 .add(simulation_options);
  
  po::variables_map vm;
  store(po::parse_command_line(argc, argv, cmdline_options), vm);
  notify(vm);
  
  __verbose      = vm.count("verbose");
  __brief_values = vm.count("brief-values");
  __print_spins  = vm.count("print-spins");
  
  /// generic options
  if ( vm.count("help") ) {
    cout << cmdline_options << endl;
    return 0;
  }
  
  if ( !L.size() ) {
    cerr << "L is required" << endl;
    return 1;
  }
  
  Assert(n_sweeps);
  n_sweeps = exp2i(lround(log2(n_sweeps)));
  if ( !vm.count("thermalization") )
    thermalization = n_sweeps/4;
  if ( !vm.count("annealing") )
    annealing = thermalization/4;
  
  Assert( JJ_IO.f == 0 || n_layer_dims > 0 );
  
  auto overwrite_ok = [&](std::ifstream &file) {
    if ( file.good() ) {
      const string pattern = "\"sweeps\" -> ";
      string line;
      bool found;
      while ( !(found = (line.compare(0, pattern.length(), pattern)==0)) && file.good() )
        getline(file, line);
      Assert_( found && file.good(), file_name );
      const size_t comma_pos = line.find(',');
      Assert( comma_pos != string::npos );
      const uint64_t old_n_sweeps = from_string<uint64_t>( line.substr(pattern.length(), comma_pos-pattern.length()) );
      if ( old_n_sweeps >= n_sweeps ) {
        cout << " already exists with " << old_n_sweeps << " sweeps" << endl;
        return false;
      }
    }
    return true;
  };
  
  if ( replace != "yes" )
  {
    std::ifstream file(file_name);
    
    if ( replace == "no" )
      Assert( !file.good() );
    else if ( replace == "auto" || replace == "auto-error" ) {
      if ( !overwrite_ok(file) )
        return replace == "auto-error";
    } else
      Assert(false);
  }
  
  vector<Float> J;
  for (IOFloat j_IO: J_IO)
    J.push_back(j_IO.f);
  
  vector<Float> J_anisotropy;
  if ( vm.count("J-anisotropy") ) {
    Assert( J_anisotropy_IO.size() == L.size() - n_layer_dims );
    FOR (d, n_layer_dims)
      J_anisotropy.push_back(0);
    for (IOFloat j_IO: J_anisotropy_IO)
      J_anisotropy.push_back(j_IO.f);
  }
  
  unique_ptr<SpinFunc>  potential;
  SpinFlipper          *spin_flipper = nullptr;
  if ( vm.count("potential") )
    Assert(false);
  else
    Assert(J.size() == 1);
  
  Assert(n);
  unique_ptr<MC> mc = MC::make(L, n);
  
  #ifndef USE_RdRand
  random_seed = std::chrono::system_clock::now().time_since_epoch().count() * getpid();
  random_engine.seed(random_seed);
  #endif
  
  Assert(J.size());
  mc->set_layer_dims(n_layer_dims);
  mc->J_ = J[0];
  if ( J_anisotropy.size() )
    mc->J_anisotropy = J_anisotropy;
  mc->JJ = JJ_IO.f;
  mc->set_update_method(update_method);
  mc->set_potential(potential.get(), potential_name);
  mc->set_flipper(spin_flipper);
  mc->set_thermalization(thermalization);
  mc->set_annealing(annealing);
  mc->set_sweeps(n_sweeps);
  
  if      ( initial_state == "low-temp" )
    mc->clear_spins();
  else if ( initial_state == "high-temp" )
    mc->randomize_spins();
  else
    Assert(false);
  
  mc->sweep();
  
  // check once more, just in case
  if ( replace == "auto" ) {
    std::ifstream file(file_name);
    if ( !overwrite_ok(file) )
      return 0;
  }
  
  if ( vm.count("file") ) {
    if ( file_name.back() == '/' ) {
      std::ostringstream ss;
      ss << file_name << " n=" << n << " L=" << L << " J=" << J << " init=" << initial_state << " updates=" << update_method;
      file_name = ss.str();
    }
    std::ofstream file(file_name);
    file << *mc;
    Assert( file.good() );
  } else
    cout << *mc << endl;
  
  if ( debug_num )
    cout << debug_num << endl;
  
  return 0;
}
