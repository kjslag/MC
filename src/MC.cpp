typedef unsigned uint;
typedef double Float;
typedef long double LongFloat;

static bool __verbose = false;
long double debug_num = 0;

//#define USE_RdRand

#include "util.hh"

typedef IO_float<Float> IOFloat;

constexpr Float pi = boost::math::constants::pi<Float>();

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
  
  Value(LongFloat min, LongFloat max, bins=1u<<7)
   : _n_samples(0), _means_k(0),
     _histo(1u<<7),
     _histo_min(min), _histo_max(max),
     _min_sample(+infinity<LongFloat>()), _max_sample(-infinity<LongFloat>())
  { _means.reserve(max_n_means); }
  
  void add(LongFloat x, LongFLoat xh)
  {
    add(x);
    
    if ( xh < _min_sample ) _min_sample = xh;
    if ( xh > _max_sample ) _max_sample = xh;
    
    int n = _histo.size();
    Assert(n);
    LongFloat bin = n*(xh-_histo_min)/(_histo_max-_histo_min);
    ++_histo[max(0,min<int>( n-1, lround(bin-.5) ))];
  }
  
  void operator<<(LongFloat x)
  { Assert(!_histo.size()); add(x) }
  
  void add(LongFloat x)
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
    return sqrt(max(var, 2*mean2)) / n_x2_samples_;
    // don't trust it, perhaps all spins are always up
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
  
  bool full(uint k) const { return _n_samples & exp2i(k); }
  
  //friend ostream& operator<<(ostream &os, const Sum &sum);
  
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

//ostream& operator<<(ostream &os, const Value::Sum &sum)
//{ return os << "sum[" << sum.x << ", " << sum.x2 << "]"; }

ostream& operator<<(ostream &os, const Value &value)
{
  os << "value["                   << value.mean()
     << ", "                       << value.error()
     << ", \"samples\" -> "        << value.n_samples();
  if ( !Value::brief ) {
  os << ",\n\"means\" -> "         << value.means()
     << ",\n\"binned errors\" -> " << value.binned_errors();
  if ( value._histo.size() ) {
  os << ",\n\"histo\" -> "         << value._histo()
     << ",\n\"histo min\" -> "     << value._histo_min
     << ",\n\"histo max\" -> "     << value._histo_max
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
  void         set_potential(const SpinFunc *V, string potential_name) {
    _potential_name = potential_name;
    set_potential_func(V);
  }
  
  virtual const SpinFunc* potential() const =0;
  
  void set_update_method(const string method);
  void set_thermalization(uint64_t thermalization) { _thermalization = thermalization; }
  void set_annealing(uint64_t annealing) { _annealing = annealing; }
  void set_sweeps(uint64_t nSweeps) { _n_sweeps = nSweeps; }
  
  void sweep()
  {
    Assert(_n_sweeps);
    
    const uint64_t nSweeps = _n_sweeps + _thermalization;
    for ( ; _sweep_num<nSweeps; ++_sweep_num ) {
      _thermalizing = _sweep_num < _thermalization;
      switch (_update_method) {
        case UpdateMethod::local:
          local_update(!_thermalizing, false);
          break;
        case UpdateMethod::global:
          if (n_ == 1)
            global_update(!_thermalizing);
          else {
            global_update(false);
            if ( potential() )
              local_update(false);
            if ( !_thermalizing )
              global_update(true);
          }
          break;
        default:
          Assert_(false, _update_method);
      }
      
      if ( (_sweep_num%16)==0 )
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
     _thermalizing(false),
     _thermalization(0),
     _annealing(0),
     _n_sweeps(0),
     _sweep_num(0),
     _histo_1 (1u<<7),
     _histo_a4(1u<<7)
  { }
  
  virtual void set_potential_func(const SpinFunc *V) =0;
  
  virtual void  local_update(bool measure, bool allow_simple=true) =0;
  virtual void global_update(bool measure) =0;
  
  UpdateMethod _update_method;
  
  std::uniform_int_distribution<Index> index_dist;
  
  string           _potential_name;
  bool             _thermalizing;
  uint64_t         _thermalization, _annealing, _n_sweeps, _sweep_num;
  Value            _sum_1, _sum_2, _sum_2_q1, _sum_2_q2, _sum_4, _sum_6, _sum_a4, _sum_cos3, _sum_tri_u, _sum_tri_v;
  vector<uint64_t> _histo_1, _histo_a4;
  Value            _n_measurement_clusters, _n_potential_clusters;
  vector<Value>    _n_clusters;
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
  
  _float norm() const { return sqrt(sq(_this)); }
  
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
_float sq(Spin_<n,_float> s) { return s|s; }

template<uint n, typename _float>
ostream& operator<<(ostream &os, const Spin_<n,_float> &s)
{ return os << s._s; }

// SpinFlipper

struct SpinFlipper {
  virtual ~SpinFlipper();
  virtual void reset() =0;
  
  virtual uint  n_flippers() const { return 1; }
  virtual uint flipper_num() const { return 0; }
  virtual bool allow_invert() const { return true; }
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

#include "SpinOperatorData.hh"

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
};
SpinFunc::~SpinFunc() {}

template<uint n>
struct SpinFunc_ : public SpinFunc {
  virtual ~SpinFunc_() {}
  virtual Float operator()(Spin_<n> s) const =0;
  virtual Spin_<n> ideal_spin() const = 0;
};

template<uint n, typename float_>
float_ s4_potential(Spin_<n,float_> s)
{ float_ V=0; FOR(a,n) V += Pow<4>(s[a]); return V; }

template<uint n>
struct s4_Potential_ : public SpinFunc_<n> {
  typedef Spin_<n> Spin;
  s4_Potential_(Float u_) : u(u_) {}
  virtual Float operator()(Spin s) const final { return u*s4_potential(s); }
  virtual vector<Float> coefficients() const final { return {u}; }
  virtual Spin ideal_spin() const final {
    Spin s;
    FOR(a,n) s[a] = u < 0 ? (a==0) : 1/sqrt(Float(n));
    return s;
  }
  const Float u;
};

unique_ptr<SpinFunc> make_s4_Potential(Float u, uint n)
{
  SpinFunc *potential = nullptr;
  if      (n==2) potential = new s4_Potential_<2>(u);
  else if (n==3) potential = new s4_Potential_<3>(u);
  else if (n==4) potential = new s4_Potential_<4>(u);
  
  Assert_(potential, n);
  return unique_ptr<SpinFunc>(potential);
}

// 2*cos(m*theta)
template<uint m, uint n, typename float_>
float_ cos_potential(Spin_<n,float_> s) {
  Assert(n==2);
  complex<Float> z(s[0], s[1]);
  return real(Pow<m>(z));
}

template<uint m>
struct cos_Potential_ : public SpinFunc_<2> {
  typedef Spin_<2> Spin;
  cos_Potential_(Float u_) : u(u_) {}
  virtual Float operator()(Spin s) const final { return u*cos_potential<m>(s); }
  virtual vector<Float> coefficients() const final { return {u}; }
  virtual Spin ideal_spin() const final
  { return u<0 ? array<Float,2>({1,0}) : array<Float,2>({cos(pi*m), sin(pi*m)}); }
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

  Assert_(potential, m);
  return unique_ptr<SpinFunc>(potential);
}

struct VisonSquare_sVBS_Potential : public SpinFunc_<4>
{
  typedef Spin_<4> Spin;
  
  VisonSquare_sVBS_Potential(Float u_, Float v_) : u(u_), v(v_) {}
  
  virtual Float operator()(Spin s) const final {
    FOR(a,4) s[a] = s[a]*s[a];
    return   u*sq(s) + v*(s[0] + s[1])*(s[2] + s[3]);
  }
  
  virtual vector<Float> coefficients() const final { return {u,v}; }
  
  virtual Spin ideal_spin() const final { Assert(false); }
  
  const Float u, v;
};

template<uint n, typename float_>
pair<float_,float_> visonTriangle_potential(Spin_<n,float_> s) {
  Assert(n==6);
  Spin_<n,float_> s2;
  FOR(a,6) s2[a] = s[a]*s[a];
  return std::make_pair(
    sq(s2[0]+s2[1])     +     sq(s2[2]+s2[3])     +     sq(s2[4]+s2[5]),
      (s2[0]-s2[1])*s[2]*s[3] + (s2[2]-s2[3])*s[4]*s[5] + (s2[4]-s2[5])*s[0]*s[1] );
}

struct VisonTrianglePotential : public SpinFunc_<6>
{
  typedef Spin_<6> Spin;
  
  VisonTrianglePotential(Float u_, Float v_) : u(u_), v(v_) {}
  
  virtual Float operator()(Spin s) const final {
    pair<Float,Float> V = visonTriangle_potential(s);
    return u*V.first + v*V.second;
  }
  
  virtual vector<Float> coefficients() const final { return {u,v}; }
  
  virtual Spin ideal_spin() const final
  {
    const Float c=cos(pi/8), s=sin(pi/8), r=sqrt(.5);
    Spin candidates[] = {
      {{1,0,0,0,0,0}}, // c-VBS
      {{c,s,0,0,0,0}}, // plaquette VBS
      {{r,r,0,0,1,0}}, // s-VBS
      {{r,r,0,0,0,1}}, // caterpillar VBS
      {{c,s,c,s,c,s}}, // star VBS
      {{s,c,s,c,s,c}}  // triangle VBS
    };
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
  
  const Float u, v;
};

// MC_

template<uint dim, // # of spacial dimensions
         uint n>   // # of spin components
class MC_ : public MC
{
  enum class InitTemp { none, high, low };
public:
  typedef Spin_<n,    Float> Spin;
  typedef Spin_<n,LongFloat> SpinSum;
  typedef array<Index,2*dim> NearestNeighbors;
  
  // POD is zero initialized by {}
  static_assert( std::is_pod<Spin>::value, "Spin isn't POD" );
  
  // lattice
  
  typedef array<uint,dim> Pos;
  
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
  
  virtual void local_update(const bool measure_, const bool allow_simple) final __attribute__((hot))
  {
    const Size n_updaes = 1+index_dist(random_engine);
    for (Size count=0; count<n_updaes; ++count) {
      const Index i      = index_dist(random_engine);
      const bool  simple = n>1 && allow_simple && bool_dist(random_engine);
      
      const Spin s1 = _spins[i];
      Spin sum_nn{};
      for(Index nn : nearestNeighbors(i))
        sum_nn += _spins[nn];
      
      Spin s2;
      Float delta_E = 0;
      if (!simple) {
        s2 = n==1 ? -s1 : Spin::random();
        delta_E += -J*((s2-s1)|sum_nn);
      } else
        s2 = (2*(s1|sum_nn)/sq(sum_nn))*sum_nn - s1;
      
      if (_V)
        delta_E += V(s2) - V(s1);
      
      delta_E = anneal(delta_E);
      if ( delta_E <= 0 || uniform_dist(random_engine) < exp(-delta_E) )
        _spins[i] = s2;
    }
    
    if ( measure_ )
      measure();
  }
  
  virtual void global_update(const bool measure_) final __attribute__((hot))
  {
    if ( _setup_global_update )
    {
      _n_clusters.resize( _flipper->n_flippers() );
      
      if ( _V ) {
        const Spin ideal_spin = _V->ideal_spin();
        const Float V_min = V(ideal_spin);
        FOR(i, 100) {
          const bool measure0 = i<3 && _flipper->allow_invert();
          _flipper->reset();
          const Spin  s1 = Spin::random();
          const Spin  s2 = _flipper->flipped(s1,measure0);
          const Spin  s3 = _flipper->flipped(s2,measure0);
          const Float V1 = V(s1);
          const Float delta_V = V(s2) - V1;
          const Float norm_s13 = (s1-s3).norm();
          Assert_( is_small(delta_V), s1, delta_V );
          Assert_( is_small(norm_s13), s1, s3, norm_s13 );
          Assert_( V_min <= V1, V_min, V1, ideal_spin, s1 );
        }
      }
      
      _setup_global_update = false;
    }
    
    if ( measure_ && !_flipper->allow_invert() ) {
      // s -> -s improved estimator code isn't allowed
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
    if (!_cluster[i0]) {
      const bool flipCluster = !force_wolff && bool_dist(random_engine);
      Index *newIndex = _newIndexStack.get();
      
      #define sum_cluster_q(op, i) do{ \
        _clusterSums[nClusters] op (s); \
        const Pos p = pos(i); \
        FOR(d, dim) \
        if ( L[d] == L[0] ) { \
          const uint x1 =    p[d]; \
          const uint x2 = (2*p[d])%L[0]; \
          _clusterSums_q[d][nClusters][0] op (_cos[x1] * s); \
          _clusterSums_q[d][nClusters][1] op (_sin[x1] * s); \
          _clusterSums_q[d][nClusters][2] op (_cos[x2] * s); \
          _clusterSums_q[d][nClusters][3] op (_sin[x2] * s); \
        } \
      } while(0)
      
      *newIndex = i0;
      _cluster[i0] = true;
      uint n_cluster_indices = 0;
      Float  cluster_delta_V = 0;
      
      {
        const Spin s = _spins[i0];
        if ( measure_ )
          sum_cluster_q(=SpinSum, i0);
        if ( force_wolff ) {
          _cluster_indices[n_cluster_indices++] = i0;
          cluster_delta_V += V(_wolff_flipper.flipped(s,false)) - V(s);
        }
        if ( flipCluster )
          flipper.flip(_spins[i0], measure_);
      }
      
      do {
        const Index j = *(newIndex--);
        
        for(Index i : nearestNeighbors(j))
          if ( !_cluster[i] ) {
            const Spin s = _spins[i];
            const Spin flipped_s = flipper.flipped(s, measure_);
            const Float delta_E = anneal( -J*(1-2*flipCluster) * ((flipped_s-s)|_spins[j]) );
            if ( delta_E > 0 && uniform_dist(random_engine) > exp(-delta_E) ) {
              *(++newIndex) = i;
              _cluster[i] = true;
              if ( measure_ )
                sum_cluster_q(+=, i);
              if ( force_wolff ) {
                _cluster_indices[n_cluster_indices++] = i;
                cluster_delta_V += V(flipped_s) - V(s);
              }
              if ( flipCluster )
                _spins[i] = flipped_s;
            }
          }
      }
      while ( newIndex+1 != _newIndexStack.get() );
      
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
      typedef SpinMatrix_<n> SpinMatrix;
      static_assert( std::is_pod<SpinMatrix>::value, "SpinMatrix isn't POD" );
      
      LongFloat  S6=0, S4_1=0;
      SpinSum    S1{}, Q2{};
      SpinMatrix M2{}, M4{};
      LongFloat  sum_2_q1=0, sum_2_q2=0;
      FOR(c, nClusters)
      {
        const SpinSum   s  = _clusterSums[c];
        const LongFloat s2 = s|s;
        
        S1 += s;
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
        if ( L[d] == L[0] )
        FOR(cs, 2) {
          sum_2_q1 += sq(_clusterSums_q[d][c][0+cs]);
          sum_2_q2 += sq(_clusterSums_q[d][c][2+cs]);
        }
      }
      const SpinMatrix M2_2 = M2|M2;
      
      uint n_q_dim = 0;
      FOR(d, dim)
        if ( L[d] == L[0] )
          ++n_q_dim;
      
      const LongFloat
      N1      = LongFloat(N),
      N2      = sq(N1),
      N4      = sq(N2),
      N6      = N2*N4,
      S1_     =    S1.norm(),
      S2      = tr(M2),
      S4      = tr(M4),
      S_2_2   = tr(M2_2),
      S_2_2_2 = tr(M2_2, M2),
      S_4_2   = tr(M4, M2),
      Q2_2    = Q2|Q2;
      
      LongFloat sum_4;
      _sum_1    << (       S1_ )/N;
      _sum_2    << (       S2  )/N2;
      _sum_2_q1 << (       sum_2_q1 )/(n_q_dim*N2);
      _sum_2_q2 << (       sum_2_q2 )/(n_q_dim*N2);
      _sum_4    << (sum_4= sq(S2) - 2*S4 + 2*S_2_2 )/N4;
      _sum_6    << (       Pow<3>(S2) - 6*S2*S4 + 16*S6 + 6*S2*S_2_2 - 24*S_4_2 + 8*S_2_2_2 )/N6;
      fill_histo(_histo_1, S1_/N);
      if (n > 1) {
        const LongFloat
        a4f = (LongFloat(3)/(2+n)),
        sum_a4_pos = (3*Q2_2 - 2*S4_1)/N4, // TODO
        sum_a4     = (3*Q2_2 - 2*S4_1 - a4f*sum_4)/N4;
        
        _sum_a4 << sum_a4;
        fill_histo(_histo_a4, sum_a4_pos/n); // TODO
      }
      if (n == 2)
        _sum_cos3 << cos_potential<3>(S1/N1);
      if (n == 6) { // TODO improved estimators
        pair<LongFloat,LongFloat> V_uv = visonTriangle_potential(S1);
        
        _sum_tri_u << (V_uv.first - .5*sum_4)/N4; // TODO I don't understand the .5
        _sum_tri_v <<  V_uv.second/N4;
      }
      
      _n_measurement_clusters << nClusters;
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
  
  Float V(Spin s) const { return (*_V)(s); }
  
protected:
  void measure()
  {
    SpinSum sum{}, sum_q1[dim][2]{}, sum_q2[dim][2]{};
    FOR(i, N) {
      const Spin s = _spins[i];
      sum += s;
      const Pos p = pos(i);
      
      FOR(d, dim)
      if ( L[d] == L[0] ) {
        const uint x1=   p[d];
        const uint x2=(2*p[d])%L[0];
        sum_q1[d][0] += _cos[x1] * s;
        sum_q1[d][1] += _sin[x1] * s;
        sum_q2[d][0] += _cos[x2] * s;
        sum_q2[d][1] += _sin[x2] * s;
      }
    }
    sum /= N;
    
    const LongFloat
    N2 = sq(LongFloat(N)),
    sum_2 = sq(sum),
    sum_4 = sq(sum_2);
    
    uint n_q_dim = 0;
    LongFloat sum_2_q1=0, sum_2_q2=0;
    FOR(d, dim)
      if ( L[d] == L[0] ) {
        ++n_q_dim;
        FOR(cs, 2) {
          sum_2_q1 += sq(sum_q1[d][cs]);
          sum_2_q2 += sq(sum_q2[d][cs]);
        }
      }
    sum_2_q1 /= n_q_dim*N2;
    sum_2_q2 /= n_q_dim*N2;
    
    const LongFloat sum_1 = sum.norm();
    
    _sum_1    <<         sum_1;
    _sum_2    <<         sum_2;
    _sum_2_q1 <<         sum_2_q1;
    _sum_2_q2 <<         sum_2_q2;
    _sum_4    <<         sum_4;
    _sum_6    <<  Pow<3>(sum_2);
    fill_histo(_histo_1, sum_1);
    if (n > 1) {
      const LongFloat
      a4f        = (LongFloat(3)/(2+n)),
      sum_a4_pos = s4_potential(sum), // TODO
      sum_a4     = sum_a4_pos - a4f*sum_4;
      
      _sum_a4 << sum_a4;
      fill_histo(_histo_a4, sum_a4_pos/n);
    }
    if (n == 2)
      _sum_cos3 << cos_potential<3>(sum);
    if (n == 6) {
      pair<LongFloat,LongFloat> V_uv = visonTriangle_potential(sum);
      
      _sum_tri_u << V_uv.first - .5*sum_4;
      _sum_tri_v << V_uv.second;
    }
  }
  
  // 0 < x < 1
  void fill_histo(vector<uint64_t> &histo, Float x) {
    ++histo[max(0,min<int>(histo.size()-1, lround(histo.size()*x-.5) ))];
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
  
  MC_(const array<uint,dim> &L__, const vector<uint> &L_v, Size N_)
    : MC(N_, L_v, n),
      L(L__),
     _spins(new Spin[N_]),
     _newIndexStack(new Index[N_]),
     _clusterSums(new SpinSum[N_]),
     _cluster_indices(new Index[N_]),
     _cos(new Float[L[0]]),
     _sin(new Float[L[0]]),
     _flipper(nullptr),
     _V(nullptr),
     _init_temp(InitTemp::none),
     _setup_global_update(true)
  {
    set_flipper(nullptr);
    
    FOR(d, dim)
    if ( L[d] == L[0] )
      _clusterSums_q[d].reset(new array<SpinSum,4>[N_]);
    
    FOR(x, L[0]) {
      _cos[x] = cos((2*pi*x)/L[0]);
      _sin[x] = sin((2*pi*x)/L[0]);
    }
  }
  
  virtual void set_potential_func(const SpinFunc *new_V) final
  {
    _V = dynamic_cast<const SpinFunc_<n>*>(new_V);
    
    if (new_V)
      Assert_(_V, n);
  }
  
  friend class MC;
  
  virtual ~MC_() {}
  
  const array<uint,dim> L; // lengths
  
private:
  const unique_ptr<Spin[]>       _spins;
  vector<bool>                   _cluster;
  const unique_ptr<Index[]>      _newIndexStack;
  const unique_ptr<SpinSum[]>    _clusterSums;
  unique_ptr<array<SpinSum,4>[]> _clusterSums_q[dim]; // 4 = #(q1, q2) * #(cos, sin)
  unique_ptr<Index[]>            _cluster_indices;
  unique_ptr<Float[]>            _cos, _sin;
  
  InvertSpin_<n>      _wolff_flipper;
  SpinFlipper_<n>    *_flipper;
  const SpinFunc_<n> *_V;
  
  InitTemp _init_temp;
  bool _setup_global_update;
};

bool __brief_values = false;
ostream& operator<<(ostream &os, const MC &mc)
{
  const SpinFunc *const potential = mc.potential();
  vector<IOFloat> J{mc.J};
  if (potential) {
    const vector<Float> coefficients = potential->coefficients();
    J.insert(J.end(), coefficients.begin(), coefficients.end());
  }
  
  const int old_prec = os.precision();
  os.precision(15);
  Value::brief = __brief_values;
  os << "{\n"
        "\"L\" -> " << mc.L_ << ",\n"
        "\"n\" -> " << mc.n_ << ",\n"
        "\"J\" -> " <<    J  << ",\n"
        "\"potential"         "\" -> \"" << mc._potential_name          << "\",\n"
        "\"update method"     "\" -> \"" << mc._update_method           << "\",\n"
        "\"sweeps"            "\" -> "   << mc._n_sweeps                <<   ",\n"
        "\"thermalization"    "\" -> "   << mc._thermalization          <<   ",\n"
        "\"annealing"         "\" -> "   << mc._annealing               <<   ",\n"
        "\"|S|"               "\" -> "   << mc._sum_1                   <<   ",\n"
        "\"histo_1"           "\" -> "   << mc._histo_1                 <<   ",\n"
        "\"S^2"               "\" -> "   << mc._sum_2                   <<   ",\n"
        "\"S^2_q1"            "\" -> "   << mc._sum_2_q1                <<   ",\n"
        "\"S^2_q2"            "\" -> "   << mc._sum_2_q2                <<   ",\n"
        "\"S^4"               "\" -> "   << mc._sum_4                   <<   ",\n"
        "\"S^6"               "\" -> "   << mc._sum_6                   <<   ",\n";
  if ( mc._sum_a4.   n_samples() )
    os<<"\"Sa^4"              "\" -> "   << mc._sum_a4                  <<   ",\n"  // should =0 if O(n) symmetry
        "\"histo_a4"          "\" -> "   << mc._histo_a4                <<   ",\n";
  if ( mc._sum_cos3. n_samples() )
    os<<"\"S cos3"            "\" -> "   << mc._sum_cos3                <<   ",\n";
  if ( mc._sum_tri_u.n_samples() )
    os<<"\"S tri u"           "\" -> "   << mc._sum_tri_u               <<   ",\n"  // should =0 if O(n) symmetry
        "\"S tri v"           "\" -> "   << mc._sum_tri_v               <<   ",\n"; // should =0 if O(n) symmetry
  Value::brief = true;
  os << "\"measurement clusters\" -> "   << mc._n_measurement_clusters  <<   ",\n"
        "\"potential clusters""\" -> "   << mc._n_potential_clusters    <<   ",\n"
        "\"clusters"          "\" -> "   << mc._n_clusters              <<    "\n";
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
  ELSE_TRY_MC(1,2)
//ELSE_TRY_MC(1,3)
//ELSE_TRY_MC(1,6)
  ELSE_TRY_MC(2,1)
//ELSE_TRY_MC(2,2)
//ELSE_TRY_MC(2,3)
  ELSE_TRY_MC(3,1)
  ELSE_TRY_MC(3,2)
  ELSE_TRY_MC(3,3)
  ELSE_TRY_MC(3,4)
  ELSE_TRY_MC(3,6)
//ELSE_TRY_MC(4,1)
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
      ("brief-values", "don't output detailed measurement errors");
  
  vector<uint> L;
  uint n = 0;
  vector<IOFloat> J_IO;
  string potential_name;
  po::options_description system_options("physics options");
  system_options.add_options()
      ("L",         po::value<vector<uint>>(&L)->multitoken(),       "lengths")
      ("n",         po::value<uint>(&n),                             "for an O(n) model")
      ("J",         po::value<vector<IOFloat>>(&J_IO)->multitoken(), "spin coupling and other potential coefficients")
      ("potential", po::value<string>(&potential_name),              "potential term: s^4; cos#; vison hexagon|square|triangle c|s-VBS");
  
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
      ("file",           po::value<string>(&file_name),                                "save file (or into directory/)")
      ("replace",        po::value<string>(&replace)->default_value("no"),             "replace file? yes, no, or auto");
  
  po::options_description cmdline_options;
  cmdline_options.add(generic_options)
                 .add(system_options)
                 .add(simulation_options);
  
  po::variables_map vm;
  store(po::parse_command_line(argc, argv, cmdline_options), vm);
  notify(vm);
  
  __verbose      = vm.count("verbose");
  __brief_values = vm.count("brief-values");
  
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
    else if ( replace == "auto" ) {
      if ( !overwrite_ok(file) )
        return 0;
    } else
      Assert(false);
  }
  
  vector<Float> J;
  for (IOFloat j_IO: J_IO)
    J.push_back(j_IO.f);
  
  unique_ptr<SpinFunc>  potential;
  SpinFlipper          *spin_flipper = nullptr;
  if ( vm.count("potential") )
  {
    auto check_n = [&n](const uint new_n) {
      if (n==0) n=new_n;
      else      Assert_(n==new_n, n, new_n); };
    
    string V_str = potential_name;
    
    if        ( V_str == "vison hexagon c-VBS" ) {
      V_str = "cos6";
    } else if ( V_str == "vison hexagon s-VBS" ) {
      check_n(3);
      V_str = "s^4";
    } else if ( V_str == "vison square c-VBS"  ) {
      V_str = "cos8";
    }
    
    if ( V_str == "s^4" ) {
      Assert_(J.size() == 1+1, J);
      potential = make_s4_Potential(J[1], n);
      if      (n==2) spin_flipper = &signed_permutation_flipper_2;
      else if (n==3) spin_flipper = &signed_permutation_flipper_3;
      else if (n==4) spin_flipper = &signed_permutation_flipper_4;
      else           Assert_(false, n);
    } else if ( V_str.substr(0,3) == "cos" ) {
      check_n(2);
      Assert_(J.size() == 1+1, J);
      const uint m = from_string<uint>(V_str.substr(3));
      potential = make_cos_Potential(J[1], m);
      if      (m==2) spin_flipper = &cos2_flipper;
      else if (m==3) spin_flipper = &cos3_flipper;
      else if (m==4) spin_flipper = &cos4_flipper;
      else if (m==6) spin_flipper = &cos6_flipper;
      else if (m==8) spin_flipper = &cos8_flipper;
      else           Assert_(false, m);
    } else if ( V_str == "vison square s-VBS" ) {
      check_n(4);
      Assert_(J.size() == 1+2, J);
      potential.reset(new VisonSquare_sVBS_Potential(J[1], J[2]));
      Assert_(false, V_str);
    } else if ( V_str == "vison triangle c-VBS" ) {
      check_n(6);
      Assert_(J.size() == 1+2, J);
      potential.reset(new VisonTrianglePotential(J[1], J[2]));
      spin_flipper = &vison_triangle_flipper;
    } else
      Assert_(false, V_str);
  }
  else
    Assert(J.size() == 1);
  
  Assert(n);
  unique_ptr<MC> mc = MC::make(L, n);
  
  #ifndef USE_RdRand
  random_seed = std::chrono::system_clock::now().time_since_epoch().count() * getpid();
  random_engine.seed(random_seed);
  #endif
  
  Assert(J.size());
  mc->J = J[0];
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
