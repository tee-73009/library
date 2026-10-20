#ifndef LIB_TEMPLATE_INOUT_HPP
#define LIB_TEMPLATE_INOUT_HPP
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
ostream &operator<<(ostream &a, const mint &b){
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
ostream &operator<<(ostream &a, const Mint &b){
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
    std::input(u...);
}
void print() { cout << "\n"; }
template <typename T, class... U, char sep = ' '>
void print(const T &t, const U &...u) {
    cout << t;
    if (sizeof...(u)) cout << sep;
    std::print(u...);
}
} // namespace std

namespace tee {
inline void debug_out() { cerr << "\n"; }
template <typename T, class... U>
void debug_out(const T &t, const U &...u) {
    cerr << t;
    if (sizeof...(u)) cerr << ", ";
    debug_out(u...);
}
#define debug(...) cerr << "L" << __LINE__ << " [" << #__VA_ARGS__ << "]: ", tee::debug_out(__VA_ARGS__)
}  // namespace tee

namespace fastio {

// 型判定用ヘルパー
template <typename T> struct is_vector : std::false_type {};
template <typename T, typename Alloc> struct is_vector<std::vector<T, Alloc>> : std::true_type {};

template <typename T> struct is_pair : std::false_type {};
template <typename T, typename U> struct is_pair<std::pair<T, U>> : std::true_type {};

// modint 型判定用ヘルパー (.val() メソッドを持つか判定)
template <typename T, typename = void>
struct is_modint : std::false_type {};

template <typename T>
struct is_modint<T, std::void_t<decltype(std::declval<T>().val())>> : std::true_type {};

template <typename T>
struct is_fastio_primitive : std::disjunction<
    std::is_same<T, char>,
    std::is_same<T, std::string>,
    std::is_integral<T>,
    std::is_same<T, __int128_t>,
    std::is_same<T, __uint128_t>
> {};

// 基本型の入力関数
template <class Char, std::enable_if_t<std::is_same_v<Char, char>, int> = 0>
inline Char in() { return getchar_unlocked(); }

template <class String, std::enable_if_t<std::is_same_v<String, std::string>, int> = 0>
inline std::string in() {
    char c; do { c = getchar_unlocked(); } while (isspace(c));
    std::string s;
    do { s.push_back(c); } while (not isspace(c = getchar_unlocked()));
    return s;
}

template <class Integer, std::enable_if_t<std::is_integral_v<Integer> or std::is_same_v<Integer, __int128_t> or std::is_same_v<Integer, __uint128_t>, int> = 0>
inline Integer in() {
    char c; do { c = getchar_unlocked(); } while (isspace(c));
    bool minus = false;
    if constexpr (std::is_signed_v<Integer> or std::is_same_v<Integer, __int128_t>) {
        if (c == '-') {
            minus = true;
            c = getchar_unlocked();
        }
    }
    using UnsignedInteger = std::conditional_t<std::is_same_v<Integer, __int128_t>, __uint128_t, std::make_unsigned_t<Integer>>;
    UnsignedInteger n = 0;
    do { n = n * 10 + (c - '0'); } while (not isspace(c = getchar_unlocked()));
    return minus ? -static_cast<Integer>(n) : static_cast<Integer>(n);
}

// 再帰的・汎用的な入力処理
template <typename T>
inline void read_single(T &t) {
    if constexpr (is_fastio_primitive<T>::value) {
        t = in<T>();
    } else if constexpr (is_modint<T>::value) {
        long long tmp = in<long long>();
        t = tmp;
    } else if constexpr (is_vector<T>::value) {
        for (auto &x : t) read_single(x);
    } else if constexpr (is_pair<T>::value) {
        read_single(t.first);
        read_single(t.second);
    } else {
        std::cin >> t;
    }
}

void input() {}
template <typename T, class... U>
void input(T &t, U &...u) {
    read_single(t);
    fastio::input(u...);
}

// 基本型の出力関数
template <class Char, std::enable_if_t<std::is_same_v<Char, char>, int> = 0>
inline void out(char c) { putchar_unlocked(c); }

template <class String, std::enable_if_t<std::is_same_v<String, std::string>, int> = 0>
inline void out(const std::string & s) { for (char c : s) putchar_unlocked(c); }

template <class Integer, std::enable_if_t<std::is_integral_v<Integer> or std::is_same_v<Integer, __int128_t> or std::is_same_v<Integer, __uint128_t>, int> = 0>
inline void out(Integer n) {
    char s[45];
    int i = 0;
    using UnsignedInteger = std::conditional_t<std::is_same_v<Integer, __int128_t>, __uint128_t, std::make_unsigned_t<Integer>>;
    UnsignedInteger u;
    if constexpr (std::is_signed_v<Integer> or std::is_same_v<Integer, __int128_t>) {
        if (n < 0) {
            putchar_unlocked('-');
            u = -static_cast<UnsignedInteger>(n);
        } else {
            u = static_cast<UnsignedInteger>(n);
        }
    } else {
        u = n;
    }
    do { s[i++] = u % 10; u /= 10; } while (u);
    while (i) putchar_unlocked(s[--i] + '0');
}

// 再帰的・汎用的な出力処理
template <typename T>
inline void write_single(const T &t) {
    if constexpr (is_fastio_primitive<T>::value) {
        out<T>(t);
    } else if constexpr (is_modint<T>::value) {
        out(t.val());
    } else if constexpr (is_vector<T>::value) {
        for (size_t i = 0; i < t.size(); ++i) {
            if (i > 0) putchar_unlocked(' ');
            write_single(t[i]);
        }
    } else if constexpr (is_pair<T>::value) {
        write_single(t.first);
        putchar_unlocked(' ');
        write_single(t.second);
    } else {
        std::cout << t;
    }
}

void print() { putchar_unlocked('\n'); }
template <typename T, class... U, char sep = ' '>
void print(const T &t, const U &...u) {
    write_single(t);
    if (sizeof...(u)) putchar_unlocked(sep);
    fastio::print(u...);
}

}  // namespace fastio

#endif