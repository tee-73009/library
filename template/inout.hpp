#pragma once
#include "template.hpp"
namespace std {
//pair_out
template <typename T, typename U>
ostream &operator<<(ostream &os, const pair<T, U> &p) {
    os << p.first << " " << p.second;
    return os;
}
//pair_in
template <typename T, typename U>
istream &operator>>(istream &is, pair<T, U> &p) {
    is >> p.first >> p.second;
    return is;
}
//vector_out
template <typename T>
ostream &operator<<(ostream &os, const vector<T> &v) {
    int s = (int)v.size();
    for (int i = 0; i < s; i++) os << (i ? " " : "") << v[i];
    return os;
}
//vector_in
template <typename T>
istream &operator>>(istream &is, vector<T> &v) {
    for (auto &x : v) is >> x;
    return is;
}
//__int128_t_in
istream &operator>>(istream &is, __int128_t &x) {
    string S;
    is >> S;
    x = 0;
    int flag = 0;
    for (auto &c : S) {
        if (c == '-') {
            flag = true;
            continue;
        }
        x *= 10;
        x += c - '0';
    }
    if (flag) x = -x;
    return is;
}
//__uint128_t_in
istream &operator>>(istream &is, __uint128_t &x) {
    string S;
    is >> S;
    x = 0;
    for (auto &c : S) {
        x *= 10;
        x += c - '0';
    }
    return is;
}
//__int128_t_out
ostream &operator<<(ostream &os, __int128_t x) {
    if (x == 0) return os << 0;
    if (x < 0) os << '-', x = -x;
    string S;
    while (x) S.push_back('0' + x % 10), x /= 10;
    reverse(begin(S), end(S));
    return os << S;
}
//__uint128_t_out
ostream &operator<<(ostream &os, __uint128_t x) {
    if (x == 0) return os << 0;
    string S;
    while (x) S.push_back('0' + x % 10), x /= 10;
    reverse(begin(S), end(S));
    return os << S;
}
//vector<vector>_out
template <typename T>
ostream &operator<<(ostream &os, const vector<vector<T>> &v)
{
    for (int i = 0; i < (int)v.size(); i++)
    {
        os << v[i] << endl;
    }
    return os;
}
//vector<vector<vector>>_out
template <typename T>
ostream &operator<<(ostream &os, const vector<vector<vector<T>>> &v)
{
    for (int i = 0; i < (int)v.size(); i++)
    {
        os << "i = " << i << endl;
        os << v[i];
    }
    return os;
}
//map_out
template <typename T, typename S>
ostream &operator<<(ostream &os, const map<T, S> &m)
{
    for (auto &[key, val] : m)
    {
        os << key << ":" << val << " ";
    }
    return os;
}
//set_out
template <typename T>
ostream &operator<<(ostream &os, const set<T> &st)
{
    auto itr = st.begin();
    for (int i = 0; i < (int)st.size(); i++)
    {
        os << *itr << (i + 1 != (int)st.size() ? " " : "");
        itr++;
    }
    return os;
}
//multiset_out
template <typename T>
ostream &operator<<(ostream &os, const multiset<T> &st)
{
    auto itr = st.begin();
    for (int i = 0; i < (int)st.size(); i++)
    {
        os << *itr << (i + 1 != (int)st.size() ? " " : "");
        itr++;
    }
    return os;
}
//queue_out
template <typename T>
ostream &operator<<(ostream &os, queue<T> q)
{
    while (q.size())
    {
        os << q.front() << " ";
        q.pop();
    }
    return os;
}
//deque_out
template <typename T>
ostream &operator<<(ostream &os, deque<T> q)
{
    while (q.size())
    {
        os << q.front() << " ";
        q.pop_front();
    }
    return os;
}
//stack_out
template <typename T>
ostream &operator<<(ostream &os, stack<T> st)
{
    while (st.size())
    {
        os << st.top() << " ";
        st.pop();
    }
    return os;
}
//priority_queue_out
template <class T, class Container, class Compare>
ostream &operator<<(ostream &os, priority_queue<T, Container, Compare> pq)
{
    while (pq.size())
    {
        os << pq.top() << " ";
        pq.pop();
    }
    return os;
}
#if __has_include(<atcoder/all>)
//998244353_in
istream &operator>>(istream &a, mint &b){
    long long tmp;
    a >> tmp;
    b = tmp;
    return a;
}
//998244353_out
ostream &operator<<(ostream &a, mint &b){
    a << b.val();
    return a;
}
//1000000007_in
istream &operator>>(istream &a, Mint &b){
    long long tmp;
    a >> tmp;
    b = tmp;
    return a;
}
//1000000007_out
ostream &operator<<(ostream &a, Mint &b){
    a << b.val();
    return a;
}
#endif
#if __has_include(<boost/multiprecision/cpp_dec_float.hpp>)
//Bint_in
istream &operator>>(istream &a, Bint &b){
    long long tmp;
    a >> tmp;
    b = tmp;
    return a;
}
//Real32_in
istream &operator>>(istream &a, Real32 &b){
    long double tmp;
    a >> tmp;
    b = tmp;
    return a;
}
//Real1024_in
istream &operator>>(istream &a, Real1024 &b){
    long double tmp;
    a >> tmp;
    b = tmp;
    return a;
}
#endif
#define ini(...) int __VA_ARGS__; input(__VA_ARGS__);
#define inl(...) long long __VA_ARGS__; input(__VA_ARGS__);
#define ins(...) string __VA_ARGS__; input(__VA_ARGS__);
#define inc(...) char __VA_ARGS__; input(__VA_ARGS__);
#define ing(name,size,n) Graph name(size); rep(i,n){inl(a,b);a--;b--;name[a].pb(b);name[b].pb(a);}
#define ing_on(name,size,n) Graph name(size); rep(i,n){inl(a,b);a--;b--;name[a].pb(b);}//有向
#define ing_cost(name,size,n) Graph_cost name(size); rep(i,n){inl(a,b,c);a--;b--;name[a].pb({b,c});name[b].pb({a,c});}//コスト付き
#define in1(s) for (int i = 0; i < (int)s.size();i++) input(s[i]);
#define in2(s, t) for (int i = 0; i < (int)s.size(); i++) input(s[i], t[i]);
#define in3(s, t, u) for (int i = 0; i < (int)s.size(); i++) input(s[i], t[i], u[i]);
#define in4(s, t, u, v) for (int i = 0; i < (int)s.size(); i++) input(s[i], t[i], u[i], v[i]);
void input() {}
template <typename T, class... U>
void input(T &t, U &...u) {
    cin >> t;
    input(u...);
}
void print() { cout << "\n"; }
template <typename T, class... U, char sep = ' '>
void print(const T &t, const U &...u) {
    cout << t;
    if (sizeof...(u)) cout << sep;
    print(u...);
}
}



namespace fastio {
template <class Char, std::enable_if_t<std::is_same_v<Char, char>, int> = 0>
inline Char in() { return getchar_unlocked(); }
template <class String, std::enable_if_t<std::is_same_v<String, std::string>, int> = 0>
inline std::string in() {
    char c; do { c = getchar_unlocked(); } while (isspace(c));
    std::string s;
    do { s.push_back(c); } while (not isspace(c = getchar_unlocked()));
    return s;
}
template <class Integer, std::enable_if_t<std::is_integral_v<Integer> and not std::is_same_v<Integer, char>, int> = 0>
inline Integer in() {
    char c; do { c = getchar_unlocked(); } while (isspace(c));
    if (std::is_signed<Integer>::value and c == '-') return -in<Integer>();
    Integer n = 0;
    do { n = n * 10 + c - '0'; } while (not isspace(c = getchar_unlocked()));
    return n;
}

template <class Char, std::enable_if_t<std::is_same_v<Char, char>, int> = 0>
inline void out(char c) { putchar_unlocked(c); }
template <class String, std::enable_if_t<std::is_same_v<String, std::string>, int> = 0>
inline void out(const std::string & s) { for (char c : s) putchar_unlocked(c); }
template <class Integer, std::enable_if_t<std::is_integral_v<Integer>, int> = 0>
inline void out(Integer n) {
    char s[20];
    int i = 0;
    if (std::is_signed<Integer>::value and n < 0) { putchar_unlocked('-'); n *= -1; }
    do { s[i ++] = n % 10; n /= 10; } while (n);
    while (i) putchar_unlocked(s[-- i] + '0');
}
#define ini(...) int __VA_ARGS__; input(__VA_ARGS__);
#define inl(...) long long __VA_ARGS__; input(__VA_ARGS__);
#define ins(...) string __VA_ARGS__; input(__VA_ARGS__);
#define inc(...) char __VA_ARGS__; input(__VA_ARGS__);
#define ing(name,size,n) Graph name(size); rep(i,n){inl(a,b);a--;b--;name[a].pb(b);name[b].pb(a);}
#define ing_on(name,size,n) Graph name(size); rep(i,n){inl(a,b);a--;b--;name[a].pb(b);}//有向
#define ing_cost(name,size,n) Graph_cost name(size); rep(i,n){inl(a,b,c);a--;b--;name[a].pb({b,c});name[b].pb({a,c});}//コスト付き
#define in1(s) for (int i = 0; i < (int)s.size();i++) input(s[i]);
#define in2(s, t) for (int i = 0; i < (int)s.size(); i++) input(s[i], t[i]);
#define in3(s, t, u) for (int i = 0; i < (int)s.size(); i++) input(s[i], t[i], u[i]);
#define in4(s, t, u, v) for (int i = 0; i < (int)s.size(); i++) input(s[i], t[i], u[i], v[i]);
void input() {}
template <typename T, class... U>
void input(T &t, U &...u) {
    t = in<T>();
    input(u...);
}
void print() { out<char>('\n'); }
template <typename T, class... U, char sep = ' '>
void print(const T &t, const U &...u) {
    out<T>(t);
    if (sizeof...(u)) out<char>(sep);
    print(u...);
}
}  // namespace fastio