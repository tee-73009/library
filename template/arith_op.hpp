#ifndef LIB_TEMPLATE_ARITH_OP_HPP
#define LIB_TEMPLATE_ARITH_OP_HPP
#include "template.hpp"
namespace tee {
namespace detail {
template <class T> struct is_pair_or_vector : false_type {};
template <class T, class U> struct is_pair_or_vector<pair<T, U>> : true_type {};
template <class T, class A> struct is_pair_or_vector<vector<T, A>> : true_type {};
template <class V> using scalar_t = enable_if_t<!is_pair_or_vector<V>::value, int>;
}
}
namespace std {
template <class T, class U>
pair<T, U> &operator+=(pair<T, U> &a, const pair<T, U> &b) {
    a.first += b.first, a.second += b.second;
    return a;
}
template <class T, class U>
pair<T, U> &operator-=(pair<T, U> &a, const pair<T, U> &b) {
    a.first -= b.first, a.second -= b.second;
    return a;
}
template <class T, class U>
pair<T, U> operator+(pair<T, U> a, const pair<T, U> &b) { return a += b; }
template <class T, class U>
pair<T, U> operator-(pair<T, U> a, const pair<T, U> &b) { return a -= b; }
template <class T, class U>
pair<T, U> operator-(const pair<T, U> &a) { return {-a.first, -a.second}; }
template <class T, class U, class V, tee::detail::scalar_t<V> = 0>
pair<T, U> operator*(const pair<T, U> &a, const V &k) { return {a.first * k, a.second * k}; }
template <class T, class U, class V, tee::detail::scalar_t<V> = 0>
pair<T, U> operator*(const V &k, const pair<T, U> &a) { return {k * a.first, k * a.second}; }
template <class T>
vector<T> &operator+=(vector<T> &a, const vector<T> &b) {
    assert(a.size() == b.size());
    for (size_t i = 0; i < a.size(); i++) a[i] += b[i];
    return a;
}
template <class T>
vector<T> &operator-=(vector<T> &a, const vector<T> &b) {
    assert(a.size() == b.size());
    for (size_t i = 0; i < a.size(); i++) a[i] -= b[i];
    return a;
}
template <class T>
vector<T> operator+(vector<T> a, const vector<T> &b) { return a += b; }
template <class T>
vector<T> operator-(vector<T> a, const vector<T> &b) { return a -= b; }
template <class T>
vector<T> operator-(vector<T> a) {
    for (auto &x : a) x = -x;
    return a;
}
template <class T, class V, tee::detail::scalar_t<V> = 0>
vector<T> operator*(vector<T> a, const V &k) {
    for (auto &x : a) x = x * k;
    return a;
}
template <class T, class V, tee::detail::scalar_t<V> = 0>
vector<T> operator*(const V &k, vector<T> a) {
    for (auto &x : a) x = k * x;
    return a;
}
}
#endif