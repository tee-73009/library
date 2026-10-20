// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/point_add_range_sum
#include "../template/template.hpp"
#include "../data_structure/segtree/dynamic_segtree.hpp"
ll op(ll a,ll b){return a+b;}
ll e(){return 0;}
signed main(){
    ll n,q;
    fastio::input(n,q);
    DynamicSegTree<ll,op,e> seg(n);
    rep(i,n){
        ll a;
        fastio::input(a);
        seg.set(i,a);
    }
    while(q--){
        ll t;
        fastio::input(t);
        if(t==0){
            ll p,x;
            fastio::input(p,x);
            seg.set(p,seg.get(p)+x);
        }else{
            ll l,r;
            fastio::input(l,r);
            fastio::print(seg.prod(l,r));
        }
    }
}
