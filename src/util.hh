// http://www.boost.org/doc/libs/1_46_0/doc/html/program_options.html
#include <boost/program_options.hpp>
#include <boost/math/constants/constants.hpp>

//#include <glib.h>

#include <fstream>
#include <iostream>
#include <sstream>
#include <chrono>
#include <string>
#include <vector>
#include <memory>
#include <random>
#include <tuple>
#include <functional>
#include <utility>
#include <complex>

#include <cstdio>
#include <cmath>
#include <cstdint>
#include <csignal>
#include <cstdint>
//#include <cassert>

#include <unistd.h> // getpid

using std::function;
using std::string;
using std::unique_ptr;
using std::array;
using std::vector;
using std::tuple;
using std::pair;
using std::ostream;
using std::istream;
using std::complex;
using std::abs; // so that cmath/abs can shadow cstdlib/abs
using std::min;
using std::max;

#ifdef DEBUG
constexpr bool debug = true;
#else
constexpr bool debug = false;
#endif

// compiler branch hints
#define likely(x)       __builtin_expect(!!(x),1)
#define unlikely(x)     __builtin_expect(!!(x),0)

// prints additional info on failed assertion
// eg: Assert_(5 < x && x < y, x, y)
#define Assert_F(x, F) \
  do { \
    if ( debug && unlikely(!(x)) ) { \
      std::cerr << "\nAssertion false: " __FILE__ ":" << __LINE__ << " : '" #x "'" << std::endl; \
      F; \
      exit(1); \
    } \
  } while(false)

#define Assert(x) Assert_F(x, )

#define Assert_(x, ...) \
        Assert_F(x, std::cerr << "{" #__VA_ARGS__ "} = " << std::make_tuple(__VA_ARGS__) << std::endl)

#define FOR(i,n) for (typename std::decay<decltype(n)>::type i=0; i<n; ++i)
#define _this (*this)

#define COUNT_ARGS(...) COUNT_ARGS_(__VA_ARGS__,9,8,7,6,5,4,3,2,1,0)
#define COUNT_ARGS_(x9,x8,x7,x6,x5,x4,x3,x2,x1,count,...) count

typedef const char *String;

template<typename T> T sq(T x) { return x*x; }
template<uint n, typename T> T Pow(T x) { return (n>=2 ? Pow<n/2>(x*x) : 1)*(n%2 ? x : 1); }

template<typename T>
ostream& operator<<(ostream &os, const vector<T> &v)
{
  const size_t n = v.size();
  os << '{';
  if (n) {
    FOR(i, n-1)
      os << v[i] << ", ";
    os << v[n-1];
  }
  return os << '}';
}

template<typename T, size_t n>
ostream& operator<<(ostream &os, const array<T,n> &a)
{
  os << '{';
  if (n) {
    FOR(i, n-1)
      os << a[i] << ", ";
    os << a[n-1];
  }
  return os << '}';
}

ostream& operator<<(ostream &os, const tuple<>&)
{ return os << "{}"; }

template<typename T, size_t i, size_t n>
struct TupleStreamHelper {
  static void print(ostream &os, const T &t)
  { os << std::get<i>(t) << ", ";
    TupleStreamHelper<T,i+1,n>::print(os, t); }
};

template<typename T, size_t n>
struct TupleStreamHelper<T, n, n> {
  static void print(ostream &os, const T &t) { os << std::get<n>(t); }
};

template<typename... Ts>
ostream& operator<<(ostream &os, const tuple<Ts...> &t)
{
  typedef tuple<Ts...> T;
  constexpr size_t n = std::tuple_size<T>::value;
  os << '{';
  TupleStreamHelper<T,0,n-1>::print(os, t);
  return os << '}';
}

template<class T>
string stringify(const T &x)
{
  std::ostringstream stream;
  if (!(stream << x))
    Assert_(false, x);
  return stream.str();
}

template<class T>
void from_string(T &x, const string &str, const bool fail_if_leftover_characters = true)
{
  std::istringstream stream(str);
  char c;
  if (!(stream >> x) || (fail_if_leftover_characters && stream.get(c)))
    Assert_(false, str);
}

template<class T>
T from_string(const string &str, const bool fail_if_leftover_characters = true)
{
  T x;
  from_string(x, str, fail_if_leftover_characters);
  return x;
}

// 2^x for integers
uint64_t exp2i(uint x) { return uint64_t(1)<<x; }

template<typename T>
bool is_small(T x) { return abs(x) < sqrt(std::numeric_limits<T>::epsilon()); }

template<typename T>
constexpr T infinity() { return std::numeric_limits<T>::infinity() }

template<typename F>
struct IO_float {
  IO_float() = default;
  IO_float(Float x) : f(x) {};
  
  F f;
};

template<typename F>
istream& operator>>(istream &is, IO_float<F> &f)
{
  F x;
  if (is >> x)
  {
    if ( is.peek() == '/' ) {
      char c;
      F y;
      if (is >> c >> y)
        f = x/y;
    }
    else if ( !is.bad() ) {
      is.clear();
      f = x;
    }
  }
  return is;
}

template<typename F>
ostream& operator<<(ostream &os, const IO_float<F> f)
{
  if ( llround(f.f) != f.f ) {
    uint den = 1u<<25;
    F    num = f.f * den;
    
    if ( llround(num) == num ) {
      do {
        num /= 2;
        den /= 2;
      } while ( llround(num)==num && den );
      num *= 2;
      den *= 2;
      Assert(den);
      return os << llround(num) << '/' << den;
    }
  }
  
  return os << f.f;
}

#ifdef USE_RdRand
class RdRand
{
public:
  typedef unsigned long long result_type; // may be unsigned short, int, or long long
  
  static constexpr result_type min() { return 0; }
  static constexpr result_type max() { return std::numeric_limits<result_type>::max(); }
  result_type operator()()
  {
    // http://gcc.gnu.org/onlinedocs/gcc/X86-Built_002din-Functions.html
    result_type rand;
    while ( !__builtin_ia32_rdrand64_step(&rand) ) {}
    return rand;
  }
};
#endif
