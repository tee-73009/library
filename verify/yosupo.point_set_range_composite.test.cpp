// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/point_set_range_composite
#include "../template/template.hpp"
#include "../data_structure/segtree.hpp"
using S=pll;
S op(S f,S g){return{g.first*f.first%mod,(g.first*f.second+g.second)%mod};}
S e(){return{1,0};}
signed main(){
    ll n,q;
    fastio::input(n,q);
    vector<S> a(n);
    rep(i,n)fastio::input(a[i]);
    Segtree<S,op,e> seg(a);
    while(q--){
        ll t;
        fastio::input(t);
        if(t==0){
            ll p;
            S f;
            fastio::input(p,f);
            seg.set(p,f);
        }else{
            ll l,r,x;
            fastio::input(l,r,x);
            S f=seg.prod(l,r);
            fastio::print((f.first*x+f.second)%mod);
        }
    }
}
