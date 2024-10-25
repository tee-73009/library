#pragma once
#include "template.hpp"
using namespace std;
#if __has_include(<atcoder/all>)
#include <atcoder/all>
using namespace atcoder;
using mint = modint998244353;
using Mint = modint1000000007;
using vm = vector<mint>;
using vvm = vector<vm>;
using vM = vector<Mint>;
using vvM = vector<vM>;
#endif
#if __has_include(<boost/multiprecision/cpp_dec_float.hpp>)
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/rational.hpp>
#include <boost/array.hpp>
namespace mp = boost::multiprecision;
// 任意長整数型
using Bint = mp::cpp_int;
// 仮数部が32桁(10進数)の浮動小数点数型
using Real32 = mp::number<mp::cpp_dec_float<32>>;
// 仮数部が1024桁(10進数)の浮動小数点数型
using Real1024 = mp::number<mp::cpp_dec_float<1024>>;
//有理数型
using Rat = boost::rational<Bint>;
#endif
#define OVERLOAD_3(_1, _2, _3, name, ...) name
#define REP1(i, n) for (auto i = std::decay_t<decltype(n)>{}; (i) < (n); ++(i))
#define REP2(i, l, r) for (auto i = (l); (i) < (r); ++(i))
#define rep(...) OVERLOAD_3(__VA_ARGS__, REP2, REP1)(__VA_ARGS__)
#define RREP1(i, n) for (auto i = (n); (i) != (std::decay_t<decltype(n)>{}); --(i))
#define RREP2(i, l, r) for (auto i = (l); (i) != (r); --(i))
#define rrep(...) OVERLOAD_3(__VA_ARGS__, RREP2, RREP1)(__VA_ARGS__)
#define rep_(i,a) rep(i,a.size())
#define all(x) std::begin(x), std::end(x)
#define rall(x) std::rbegin(x), std::rend(x)
#define yesno(T) if(T){cout<<"yes"<<endl;}else{cout<<"no"<<endl;}
#define YesNo(T) if(T){cout<<"Yes"<<endl;}else{cout<<"No"<<endl;}
#define YESNO(T) if(T){cout<<"YES"<<endl;}else{cout<<"NO"<<endl;}
#define fin(...) ({print(__VA_ARGS__); return 0; })
#define break_with(...) ({ __VA_ARGS__; break; })
#define continue_with(...) ({ __VA_ARGS__; continue; })
using ll=long long;
using ld=long double;
using ull=unsigned long long;
using pii = pair<int, int>;
using pll = pair<ll, ll>;
using vi = vector<int>;
using vll = vector<long long>;
using vs = vector<string>;
using vc = vector<char>;
using vb = vector<bool>;
using vvi = vector<vi>;
using vvll = vector<vll>;
using vvs = vector<vs>;
using vvc = vector<vc>;
using vvb = vector<vb>;
using vvvi = vector<vvi>;
using vvvll = vector<vvll>;
using vvvb = vector<vvb>;
using vpll = vector<pll>;
using mi = map<int, int>;
using mll = map<ll, ll>;
using msi = map<string, int>;
using msl = map<string, ll>;
using mmi = map<int,mi>;
using mml = map<ll,mll>;
using Graph = vector<vector<ll>>;
const string ABC="ABCDEFGHIJKLMNOPQRSTUVWXYZ";
const string abc="abcdefghijklmnopqrstuvwxyz";
#define gcd __gcd
#define pb push_back
template<class T> using pq = priority_queue<T>;
template<class T> using pqr = priority_queue<T, vector<T>, greater<T>>;
#define lb(v, k) (lower_bound((v).begin(), (v).end(), (k)) - v.begin())
#define ub(v, k) (upper_bound((v).begin(), (v).end(), (k)) - v.begin())
#define fi first
#define se second
#define elif else if
#define updiv(n,x) (n + x - 1) / x
#define rounddiv(n,x) (ll)((double)(n)/(double)(x)+0.5)
#define fix(n) fixed << setprecision(n)
#define _print(n) cout << " " << (n) << endl
#define print_(n) cout << (n) << " "
#define ioinit ios::sync_with_stdio(false);std::cin.tie(nullptr)
template<class... T>
//3つ以上でも使えるmax,min
constexpr auto min(T... a){
    return min(initializer_list<common_type_t<T...>>{a...});
}
template<class... T>
constexpr auto max(T... a){
    return max(initializer_list<common_type_t<T...>>{a...});
}
template <typename T>
inline bool chmax(T &a, T b) { return ((a < b) ? (a = b, true) : (false)); }
template <typename T>
inline bool chmin(T &a, T b) { return ((a > b) ? (a = b, true) : (false)); }
const ll mod=998244353;
const ll MOD=1000000007;
const ld PI=3.141592653589793;
const vll dx4={0,-1,0,1};
const vll dy4={1,0,-1,0};
const vll dx8={1,0,-1,0,1,-1,1,-1};
const vll dy8={0,1,0,-1,1,-1,-1,1};
const vll dx6={1,0,-1,0,1,-1};
const vll dy6={0,1,0,-1,1,-1};
const int inf = INT_MAX / 2;
const ll INF= 1LL << 60;
