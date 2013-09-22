// make && echo > jobs && for L in 64 128; do for beta in 0 0.25 0.5 1 2; do f="results/ising_${L}_${beta}"; echo -e "$f:\n\t./MC --L $L --L $L --beta $beta --sweep 10000 --file $f\n" >> jobs; done; done

typedef unsigned uint;
typedef double Float;

static bool verbose = false;
double debug_num = 0;

#include "util.hh"

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
  
  Float              beta;
  const Size         N; // number of spins
  const vector<uint> L_;
  const uint         n_;
  
  friend std::ostream& operator<<(std::ostream &os, const MC &mc);
  
protected:
  MC(Size N_, const vector<uint> &L_v, uint n__)
    : beta(0),
      N(N_),
      L_(L_v),
      n_(n__),
     _newSpinStack(new Index[N_]),
      index_dist(0, N_-1),
      _nSweeps(0)
  { }
  
  vector<bool>              _cluster;
  const unique_ptr<Index[]> _newSpinStack;
  Index *_newSpin;
  Size _nFlip;
  
  std::mt19937                          random_engine; // todo seed
  std::uniform_int_distribution<Index>  index_dist;
  std::uniform_real_distribution<Float> uniform_dist;
  std::normal_distribution<Float>       normal_dist;
  
  uint64_t _nSweeps;
  vector<double> _sum1, _sum2, _sum4;
};
MC::~MC() {}

ostream& operator<<(ostream &os, const MC &mc)
{
  os << "{\n"
     << "\"n\" -> "    << mc.n_ << ",\n"
     << "\"L\" -> "    << mc.L_ << ",\n"
     << "\"beta\" -> " << mc.beta << ",\n"
     << "\"moments\" -> {" << mc._nSweeps << ", "
                           << mc._sum1    << ", "
                           << mc._sum2    << ", "
                           << mc._sum4    << "}\n";
  return os << "}\n";
}

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
    
    ++_nSweeps;
    _sum1.push_back( sqrt(avg2) );
    _sum2.push_back( avg2       );
    _sum4.push_back( avg2*avg2  );
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
  
  MC_(const array<uint,dim> &L_, const vector<uint> &L_v, const array<uint,dim> &Lp, Size N_)
    : MC(N_, L_v, n),
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
    mc = new MC_<dim_,n_>(L_, L, Lp, N); \
  }
  
  if (false) {}
  ELSE_TRY_MC(1,1)
  //ELSE_TRY_MC(1,2)
  ELSE_TRY_MC(2,1)
  //ELSE_TRY_MC(2,2)
  //ELSE_TRY_MC(3,6)
  else
    Assert(false, L, N);
  
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
  namespace po = boost::program_options;
  
  po::options_description generic_options("generic options");
  generic_options.add_options()
      ("help,h",                               "print help message")
      ("verbose,v",                            "print verbose output");
  
  po::options_description system_options("physics options");
  system_options.add_options()
      ("L",    po::value<vector<uint>>(),             "lengths")
      ("n",    po::value<uint>()->default_value(1),   "for an O(n) model")
      ("beta", po::value<double>()->default_value(0), "thermodynamic beta");
  
  po::options_description simulation_options("simulation options");
  simulation_options.add_options()
      ("sweeps", po::value<uint64_t>()->default_value(1), "# of MC sweeps")
      ("file",  po::value<string>(),                      "save file");
  
  po::options_description cmdline_options;
  cmdline_options.add(generic_options)
                 .add(system_options)
                 .add(simulation_options);
  
  po::variables_map vm;
  store(po::parse_command_line(argc, argv, cmdline_options), vm);
  notify(vm);
  
  verbose = vm.count("verbose"); // todo
  
  /// generic options
  if ( vm.count("help") )
  {
    std::cout << cmdline_options << std::endl;
    return 0;
  }
  
  if ( !vm.count("L") ) // todo: file
  {
    std::cerr << "L is required" << std::endl;
    return 1;
  }
  
  MC *mc = MC::make(vm["L"].as<vector<uint>>(), vm["n"].as<uint>());
  mc->clear_spins();
  mc->beta = vm["beta"].as<double>();
  
  const uint64_t nSweeps = vm["sweeps"].as<uint64_t>();
  for (uint64_t sweep=0; sweep<nSweeps; ++sweep)
    mc->update();
  
  if ( vm.count("file") )
  {
    std::ofstream file( vm["file"].as<string>() );
    file << *mc;
  }
  
  delete mc;
  
  std::cout << debug_num << std::endl;
  
  return 0;
}
