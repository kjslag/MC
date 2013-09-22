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
  do \
  { \
    if ( debug && unlikely(!(x)) ) \
    { \
      print_backtrace(); \
      std::cerr << "\nAssertion false: " __FILE__ ":" << __LINE__ << " : '"#x"'" << std::endl \
                << #__VA_ARGS__ " = " << std::make_tuple(__VA_ARGS__) << std::endl; \
      exit(1); \
    } \
  } while(0)

#define FOR(i,n) for (uint i=0; i<n; ++i)
#define _this (*this)

/*
template<size_t i, size_t n, typename... Ts>
ostream& tuple_streamer_helper(ostream &os, const tuple<Ts...> &t)
{ return os << std::get<n-i>(t) << ", "; }

template<typename... Ts>
ostream& tuple_streamer_helper<0,0,Ts...>(ostream &os, const tuple<Ts...> &)
{ return os; }

template<size_t n, typename... Ts>
ostream& tuple_streamer_helper<1,n,Ts...>(ostream &os, const tuple<Ts...> &t)
{ return os << std::get<n-1>(t); }

template<typename... Ts>
ostream& operator<<(ostream &os, const tuple<Ts> &t)
{
  constexpr n = std::tuple_size<decltype(t)>::value;
  return os << '{' << tuple_streamer_helper<n,n>(os, t) << '}';
}
*/

template<typename T>
ostream& operator<<(ostream &os, const vector<T> &v)
{
  const size_t n = v.size();
  os << '{';
  FOR(i, n-1)
    os << v[i] << ", ";
  return os << v[n-1] << '}';
}

template<typename T, size_t n>
ostream& operator<<(ostream &os, const array<T,n> &a)
{
  os << '{';
  FOR(i, n-1)
    os << a[i] << ", ";
  return os << a[n-1] << '}';
}

ostream& operator<<(ostream &os, const tuple<>&)
{ return os << "{}"; }

template<typename T0>
ostream& operator<<(ostream &os, const tuple<T0> &t)
{ return os << '{' << std::get<0>(t) << '}'; }

template<typename T0, typename T1>
ostream& operator<<(ostream &os, const tuple<T0,T1> &t)
{ return os << '{' << std::get<0>(t) << ", " << std::get<1>(t) << '}'; }

template<typename T0, typename T1, typename T2>
ostream& operator<<(ostream &os, const tuple<T0,T1,T2> &t)
{ return os << '{' << std::get<0>(t) << ", " << std::get<1>(t) << ", " << std::get<2>(t) << '}'; }

template<typename T0, typename T1, typename T2, typename T3>
ostream& operator<<(ostream &os, const tuple<T0,T1,T2,T3> &t)
{ return os << '{' << std::get<0>(t) << ", " << std::get<1>(t) << ", " << std::get<2>(t) << ", " << std::get<3>(t) << '}'; }

template<typename T0, typename T1, typename T2, typename T3, typename T4>
ostream& operator<<(ostream &os, const tuple<T0,T1,T2,T3,T4> &t)
{ return os << '{' << std::get<0>(t) << ", " << std::get<1>(t) << ", " << std::get<2>(t) << ", " << std::get<3>(t) << ", " << std::get<4>(t) << '}'; }

template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5>
ostream& operator<<(ostream &os, const tuple<T0,T1,T2,T3,T4,T5> &t)
{ return os << '{' << std::get<0>(t) << ", " << std::get<1>(t) << ", " << std::get<2>(t) << ", " << std::get<3>(t) << ", " << std::get<4>(t) << ", " << std::get<5>(t) << '}'; }

template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
ostream& operator<<(ostream &os, const tuple<T0,T1,T2,T3,T4,T5,T6> &t)
{ return os << '{' << std::get<0>(t) << ", " << std::get<1>(t) << ", " << std::get<2>(t) << ", " << std::get<3>(t) << ", " << std::get<4>(t) << ", " << std::get<5>(t) << ", " << std::get<6>(t) << '}'; }

template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7>
ostream& operator<<(ostream &os, const tuple<T0,T1,T2,T3,T4,T5,T6,T7> &t)
{ return os << '{' << std::get<0>(t) << ", " << std::get<1>(t) << ", " << std::get<2>(t) << ", " << std::get<3>(t) << ", " << std::get<4>(t) << ", " << std::get<5>(t) << ", " << std::get<6>(t) << ", " << std::get<7>(t) << '}'; }

static void print_backtrace()
{
  //g_on_error_query("MC");
  constexpr int size = 20;
  void *buffer[size];
  backtrace_symbols_fd(buffer, backtrace(buffer, size), STDOUT_FILENO);
}
