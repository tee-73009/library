// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/range_affine_point_get
#include "../template/template.hpp"
#include "../data_structure/dual_segtree.hpp"
using S=pair<ll,ll>;
using F=pair<ll,ll>;
S op(S a,S b){return{(a.first+b.first)%mod,a.second+b.second};}
S e(){return{0,0};}
S mapping(F f,S x){return{(f.first*x.first+f.second*x.second)%mod,x.second};}
F composition(F f,F g){return{f.first*g.first%mod,(f.first*g.second+f.second)%mod};}
F id(){return{1,0};}
signed main(){
    ll n,q;
    fastio::input(n,q);
    DualSegtree<S,F,mapping,composition,id> seg(n);
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
