// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/point_add_range_sum
#include "../template/template.hpp"
#include "../data_structure/fenwicktree.hpp"
signed main() {
    ll n,q;
    fastio::input(n,q);
    BIT<ll> bit(n);
    while(q--){
        ll t;
        fastio::input(t);
        if(t==0){
            ll p,x;
            fastio::input(p,x);
            bit.add(p,x);
        }
        else{
            ll l,r;
            fastio::input(l,r);
            fastio::print(bit.sum(l,r));
        }
    }
}