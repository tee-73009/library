#define PROBLEM "https://judge.yosupo.jp/problem/dynamic_sequence_range_affine_range_sum"
#include "../template/template.hpp"
#include "../data_structure/implicit_treap.hpp"
struct S{
    mint value;
    ll size;
};
using F = pair<ll, ll>;
S op(S a, S b) { return {a.value + b.value, a.size + b.size}; }
S e() { return S{0, 0}; }
S mapping(F f, S x) { return {(f.fi * x.value + f.se * x.size), x.size}; }
F composition(F f, F g) { return {f.fi * g.fi, (f.fi * g.se + f.se)}; }
F id() { return F{1, 0}; }
int main() {
    ll n, q;
    fastio::input(n, q);
    vector<S> a(n);
    rep(i, n){
        ll x;
        fastio::input(x);
        a[i] = S{x, 1};
    }
    ImplicitTreap<S, op, e, F, mapping, composition, id> t(a);
    while (q--) {
        int ty;
        cin >> ty;
        if (ty == 0) { int i; ll x; fastio::input(i, x); t.insert(i, {x, 1}); }
        else if (ty == 1) { int i; fastio::input(i); t.erase(i); }
        else if (ty == 2) { int l, r; fastio::input(l, r); t.reverse(l, r); }
        else if (ty == 3) { int l, r; ll b, c; fastio::input(l, r, b, c); t.apply(l, r, {b, c}); }
        else { int l, r; cin >> l >> r; fastio::print(t.prod(l, r).value); }
    }
}
