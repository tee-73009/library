// competitive-verifier: PROBLEM https://onlinejudge.u-aizu.ac.jp/courses/library/3/DSL/2/DSL_2_G
#include "../template/template.hpp"
#include "../data_structure/dynamic_lazysegtree.hpp"
using S=pair<ll,ll>;  // (和,長さ)
S op(S a,S b){return{a.first+b.first,a.second+b.second};}
S e(){return{0,0};}
S mapping(ll f,S x){return{x.first+f*x.second,x.second};}
ll composition(ll f,ll g){return f+g;}
ll id(){return 0;}
signed main(){
    ll n,q;
    fastio::input(n,q);
    DynamicLazySegtree<S,op,e,ll,mapping,composition,id> seg(n,S{0,1});
    while(q--){
        ll t,s,tt;
        fastio::input(t,s,tt);
        s--;  // 1-indexed の閉区間 -> 半開区間 [s-1,tt)
        if(t==0){
            ll x;
            fastio::input(x);
            seg.apply(s,tt,x);
        }else fastio::print(seg.prod(s,tt).first);
    }
}
