typedef unsigned uint;
typedef double Float;
typedef long double LongFloat;

static bool verbose = false;
long double debug_num = 0;

//#define USE_RdRand

#include "util.hh"
#include "Value.hh"
#include "MC.hh"

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
