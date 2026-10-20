// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/point_set_range_composite_large_array
#include "../template/template.hpp"
#include "../data_structure/dynamic_segtree.hpp"
using S=pair<ll,ll>;
S op(S a,S b){return{a.first*b.first%mod,(a.second*b.first+b.second)%mod};}
S e(){return{1,0};}
signed main(){
    ll n,q;
    fastio::input(n,q);
    DynamicSegtree<S,op,e> seg(n);
    while(q--){
        ll t;
        fastio::input(t);
        if(t==0){
            ll p,c,d;
            fastio::input(p,c,d);
            seg.set(p,{c,d});
        }
        else{
            ll l,r,x;
            fastio::input(l,r,x);
            auto[c,d]=seg.prod(l,r);
            fastio::print(((c*x)%mod+d)%mod);
        }
    }
}
