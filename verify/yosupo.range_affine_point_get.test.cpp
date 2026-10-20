// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/range_affine_point_get
#include "../template/template.hpp"
#include "../data_structure/dual_segtree.hpp"
using S=ll;
using F=pair<ll,ll>;
S mapping(F f,S x){return(f.first*x%mod+f.second)%mod;}
F composition(F f,F g){return{g.first*f.first%mod,(g.first*f.second+g.second)%mod};}
F id(){return{1,0};}
signed main(){
    ll n,q;
    fastio::input(n,q);
    vll a(n);
    fastio::input(a);
    DualSegtree<S,F,mapping,composition,id> seg(a);
    rep(i,n) seg.set(i,a[i]);
    while(q--){
        ll t;
        fastio::input(t);
        if(t==0){
            ll l,r,b,c;
            fastio::input(l,r,b,c);
            seg.apply(l,r,{b,c});
        }
        else{
            ll i;
            fastio::input(i);
            fastio::print(seg.get(i));
        }
    }
}
