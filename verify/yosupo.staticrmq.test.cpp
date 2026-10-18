#define PROBLEM "https://judge.yosupo.jp/problem/staticrmq"
#include "../template/template.hpp"
#include "../data_structure/sparsetable.hpp"
ll op(ll i,ll j){
    return min(i,j);
}
signed main() {
    ll n,q;
    fastio::input(n,q);
    vll a(n);
    fastio::input(a);
    SparseTable<ll,op> table(a);
    while(q--){
        ll l,r;
        fastio::input(l,r);
        fastio::print(table.prod(l,r));
    }
}