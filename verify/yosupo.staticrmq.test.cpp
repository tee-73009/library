// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/staticrmq
#include "../template/template.hpp"
#include "../data_structure/sparsetable.hpp"
int op(int a,int b){return min(a,b);}
signed main(){
    ll n,q;
    fastio::input(n,q);
    vector<int> a(n);
    rep(i,n)fastio::input(a[i]);
    SparseTable<int,op> table(a);
    while(q--){
        ll l,r;
        fastio::input(l,r);
        fastio::print(table.prod(l,r));
    }
}
