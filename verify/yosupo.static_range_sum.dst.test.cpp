// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/static_range_sum
#include "../template/template.hpp"
#include "../data_structure/disjointsparsetable.hpp"
ll op(ll a,ll b){return a+b;}
signed main(){
    ll n,q;
    fastio::input(n,q);
    vector<ll> a(n);
    rep(i,n)fastio::input(a[i]);
    DisjointSparseTable<ll,op> table(a);
    while(q--){
        ll l,r;
        fastio::input(l,r);
        fastio::print(table.prod(l,r));
    }
}
