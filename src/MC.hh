
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
