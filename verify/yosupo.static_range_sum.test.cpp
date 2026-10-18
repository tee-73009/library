#define PROBLEM "https://judge.yosupo.jp/problem/static_range_sum"
#include "../template/template.hpp"
#include "../data_structure/disjointsparsetable.hpp"
ll op(ll a, ll b) { return a + b; }
signed main() {
    ioinit;
    ll n,q;
    fastio::input(n,q);
    vll a(n);
    fastio::input(a);
    DisjointSparseTable<ll, op> table(a);
    rep(i, q) {
        ll l,r;
        fastio::input(l,r);
        fastio::print(table.prod(l,r));
    }
}