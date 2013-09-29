// http://www.boost.org/doc/libs/1_46_0/doc/html/program_options.html
#include <boost/program_options.hpp>

//#include <glib.h>

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <random>
#include <tuple>
#include <functional>
#include <utility>

#include <cstdio>
#include <cmath>
#include <cstdint>
#include <csignal>
#include <cstdint>
//#include <cassert>

#include <execinfo.h>

using std::string;
using std::unique_ptr;
using std::array;
using std::vector;
using std::tuple;
using std::pair;
using std::ostream;
using std::istream;

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
#define Assert(x, ...) \
  do { \
    if ( debug && unlikely(!(x)) ) { \
      print_backtrace(); \
      std::cerr << "\nAssertion false: " __FILE__ ":" << __LINE__ << " : '"#x"'" << std::endl \
                << #__VA_ARGS__ " = " << std::make_tuple(__VA_ARGS__) << std::endl; \
      exit(1); \
    } \
  } while(false)

#define FOR(i,n) for (typename std::remove_const<decltype(n)>::type i=0; i<n; ++i)
#define _this (*this)

template<typename T>
ostream& operator<<(ostream &os, const vector<T> &v)
{
  const size_t n = v.size();
  os << '{';
  if (n)
  {
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
  if (n)
  {
    FOR(i, n-1)
      os << a[i] << ", ";
    os << a[n-1];
  }
  return os << '}';
}

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
  if (n)
    TupleStreamHelper<T,0,n-1>::print(os, t);
  return os << '}';
}

static void print_backtrace()
{
  //g_on_error_query("MC");
  constexpr int size = 20;
  void *buffer[size];
  backtrace_symbols_fd(buffer, backtrace(buffer, size), STDOUT_FILENO);
}
