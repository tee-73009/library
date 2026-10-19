// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/dynamic_sequence_range_affine_range_sum
#include "../template/template.hpp"
#include "../data_structure/implicit_treap.hpp"

struct S {
    mint value;
    ll size;
};
using F = pair<mint, mint>;

S op(S a, S b) { return {a.value + b.value, a.size + b.size}; }
S e() { return {0, 0}; }
S mapping(F f, S x) { return {f.first * x.value + f.second * x.size, x.size}; }
F composition(F f, F g) { return {f.first * g.first, f.first * g.second + f.second}; }
F id() { return {1, 0}; }

int main() {
    ll n, q;
    fastio::input(n, q);
    vector<S> a(n);
    rep(i, n) {
        mint x;
        fastio::input(x);
        a[i] = S{x, 1};
    }
    ImplicitTreap<S, op, e, F, mapping, composition, id> t(a);
    while (q--) {
        int ty;
        fastio::input(ty);
        if (ty == 0) {
            int i;
            mint x;
            fastio::input(i, x);
            t.insert(i, S{x, 1});
        } else if (ty == 1) {
            int i;
            fastio::input(i);
            t.erase(i);
        } else if (ty == 2) {
            int l, r;
            fastio::input(l, r);
            t.reverse(l, r);
        } else if (ty == 3) {
            int l, r;
            mint b, c;
            fastio::input(l, r, b, c);
            t.apply(l, r, {b, c});
        } else {
            int l, r;
            fastio::input(l, r);
            fastio::print(t.prod(l, r).value);
        }
    }
}