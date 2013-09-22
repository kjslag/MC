typedef unsigned uint;
typedef double Float;

double debug_num = 0;

#include "util.hh"

struct Value
{
  Value() : _count(0) {}
  
  void operator<<(double x)
  {
    uint k = 0;
    next_bin:
    
    if ( k < _bins.size() )
    {
      Bin &bin = _bins[k];
      
      bin.x2 += x*x;
      if ( !full(k) )
        bin.x = x;
      else
      {
        x += bin.x;
        ++k;
        goto next_bin;
      }
    }
    else
      _bins.push_back(x);
    
    ++_count;
  }
  
  double sum() const
  {
    double sum = 0;
    const uint n = _bins.size();
    FOR(k,n)
      if ( full(k) )
        sum += _bins[k].x;
    return sum;
  }
  
private:
  struct Bin
  {
    Bin(double y) : x(y), x2(y*y) {}
    
    double x, x2;
  };
  
  bool full(uint k) const
  { return _count & (uint64_t(1) << k); }
  
  vector<Bin> _bins;
  uint64_t    _count;
};

template<uint n>
struct Spin_
{
  Float& operator[](uint k)       { return _s[k]; }
  Float  operator[](uint k) const { return _s[k]; }
  
  void  operator+=(const Spin_ &s) { FOR(k,n) _s[k] += s[k]; }
  void  operator*=(Float        x) { FOR(k,n) _s[k] *= x;    }
  
  Float operator|(const Spin_ &s) const { Float x=0; FOR(k,n) x += _s[k]*s[k]; return x; } // dot product
  
  void flip(Spin_ r)
  { r *= Float(-2) * (_this|r);
    _this += r; }
  
  void normalize()
  { Float norm = 1/sqrt(_this | _this);
    FOR(k,n) _s[k] *= norm; }
  
private:
  array<Float,n> _s;
};

class MC
{
public:
  // lattice
  
  typedef uint Size;
  
  typedef Size Index;
  
  //
  
  static MC* make(const vector<uint> &L, uint n_fields);
  
  virtual ~MC();
  
  virtual void clear_spins() =0;
  virtual void randomize_spins() =0;
  virtual void update() =0;
  
  Float      beta;
  const Size N; // number of spins
  
protected:
  MC(Size N_)
    : beta(0),
      N(N_),
     _newSpinStack(new Index[N_]),
      index_dist(0, N_-1)
  { }
  
  vector<bool>              _cluster;
  const unique_ptr<Index[]> _newSpinStack;
  Index *_newSpin;
  Size _nFlip;
  
  std::mt19937                          random_engine; // todo seed
  std::uniform_int_distribution<Index>  index_dist;
  std::uniform_real_distribution<Float> uniform_dist;
  std::normal_distribution<Float>       normal_dist;
  
  Value _sum1, _sum2, _sum4;
};
MC::~MC() {}

std::ostream& operator<<(std::ostream &os, const MC &mc)
{ return os << "{}"; } // TODO

template<uint dim,              // # of spacial dimensions
         uint n>                // # of fields
class MC_ : public MC
{
public:
  typedef Spin_<n> Spin;
  
  // lattice
  
  typedef array<uint,dim> Pos;
  
  Index index(const Pos p) const
  {
    Index i = p[0];
    for (uint d=1; d<dim; ++d)
      i = i*L[d] + p[d];
    return i;
  }
  
  Pos pos(Index i) const
  {
    Pos p;
    for (int d=dim-1; d>=0; --d)
    {
      p[d] = i%L[d];
      i /= L[d];
    }
    return p;
  }
  
  // MC
  
  virtual ~MC_() {}
  
  virtual void clear_spins()
  {
    Spin s;
    FOR(k,n) s[k] = (k==0);
    FOR(i,N) _spins[i] = s;
  }
  
  virtual void randomize_spins()
  { FOR(i,N) _spins[i] = random_spin(); }
  
  virtual void update() __attribute__((hot))
  {
    _nFlip = 0;
    while (2*_nFlip < N)
    {
      _cluster.assign(N, false);
      
      const Spin r = random_spin();
      _newSpin = _newSpinStack.get()-1;
      add_spin(index_dist(random_engine), r);
      
      do
      {
        const Index j = *(_newSpin--);
        const Pos   q = pos(j);
        
        FOR(d, dim)
        for (int dir=-1; dir<=+1; dir+=2)
        {
          Pos p = q;
          p[d] = (p[d] + L[d] + dir) % L[d];
          const Index i = index(p);
          if ( !_cluster[i] )
          {
            const Float delta_E = Float(2)*beta*(r|_spins[i])*(r|_spins[j]);
            if ( uniform_dist(random_engine) > exp(delta_E) )
              add_spin(i, r);
          }
        }
      }
      while ( _newSpin+1 != _newSpinStack.get() );
    }
    
    array<double,n> avg = {};
    FOR(i,N) FOR(k,n) avg[k] += _spins[i][k];
    
    double avg2 = 0;
    FOR(k,n) avg2 += avg[k]*avg[k];
    
    _sum1 << sqrt(avg2);
    _sum2 << avg2;
    _sum4 << avg2*avg2;
  }
  
protected:
  void add_spin(const Index i, const Spin &r)
  {
    ++_nFlip;
    _spins[i].flip(r);
    _cluster[i] = true;
    *(++_newSpin) = i;
  }
  
  Spin random_spin()
  {
    Spin s;
    FOR(k,n) s[k] = normal_dist(random_engine);
    s.normalize();
    return s;
  }
  
  MC_(const array<uint,dim> &L_, const array<uint,dim> &Lp, Size N_)
    : MC(N_),
      L(L_),
     _Lp(Lp),
     _spins(new Spin[N_])
  { }
  
  const array<uint,dim> L; // lengths
  
  friend class MC;
  
private:
  const array<uint,dim> _Lp; // _Lp[d] = L[d+1] * ... * L[dim]
  
  const unique_ptr<Spin[]> _spins;
};

MC* MC::make(const vector<uint> &L, uint n_fields)
{
  size_t N  = 1;
  long double Nf = 1;
  const uint dim = L.size();
  FOR(d,dim) { N *= L[d]; Nf *= L[d]; }
  
  const size_t Nmax = std::numeric_limits<MC::Size>::max();
  Assert( N <= Nmax && Nf <= Nmax, N );
  
  MC *mc = nullptr;
  #define ELSE_TRY_MC(dim_, n_) \
  else if (dim == dim_ && n_fields == n_) \
  { \
    array<uint,dim_> L_, Lp; \
    uint Lp0 = 1; \
    for (int d=dim-1; d>=0; --d) { L_[d] = L[d]; Lp[d] = Lp0; Lp0 *= L[d]; } \
    mc = new MC_<dim_,n_>(L_, Lp, N); \
  }
  
  if (false) {}
  //ELSE_TRY_MC(1,1)
  //ELSE_TRY_MC(1,2)
  //ELSE_TRY_MC(2,1)
  //ELSE_TRY_MC(2,2)
  ELSE_TRY_MC(3,6)
  else
    Assert(false, dim, N);
  
  #undef ELSE_TRY_MC
  
  return mc;
}

static void single_handler(int sig)
{
  printf("Error: signal %d:\n", sig);
  
  print_backtrace();
  
  exit(1);
}

int main(const int argc, char *argv[])
{
  // http://www.gnu.org/s/hello/manual/libc/Signal-Handling.html
  signal( SIGSEGV, single_handler );
  //signal( SIGABRT, single_handler );
  
  // read program options
  if (false)
  {
    namespace po = boost::program_options;
    
    po::options_description generic_options("generic options");
    generic_options.add_options()
        ("help,h",                                              "print help message")
        ("verbose,v",                                           "print verbose output");
    
    po::options_description cmdline_options;
    cmdline_options.add(generic_options);
    
    po::variables_map vm;
    store(po::command_line_parser(argc, argv).
          options(cmdline_options).run(), vm);
    notify(vm);
  }
  
  constexpr uint L = 150;
  MC *mc = MC::make({L,L,L}, 6);
  mc->clear_spins();
  mc->beta = 1;
  mc->update();
  delete mc;
  
  std::cout << debug_num << std::endl;
  
  return 0;
}
