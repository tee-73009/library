// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/range_chmin_chmax_add_range_sum
#include "../template/template.hpp"
#include "../data_structure/segtree/segtree_beats.hpp"
signed main(){
    ll n,q;
    fastio::input(n,q);
    vector<ll> a(n);
    rep(i,n)fastio::input(a[i]);
    SegTreeBeats seg(a);
    while(q--){
        ll t,l,r;
        fastio::input(t,l,r);
        if(t==3){
            fastio::print(seg.sum(l,r));
            continue;
        }
        ll b;
        fastio::input(b);
        if(t==0)seg.chmin(l,r,b);
        else if(t==1)seg.chmax(l,r,b);
        else seg.add(l,r,b);
    }
}
