// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/static_range_sum
#include "../template/template.hpp"
#include "../util/get_sum.hpp"
signed main(){
    ll n,q;
    fastio::input(n,q);
    vector<ll> a(n);
    rep(i,n)fastio::input(a[i]);
    vector<ll> ar;
    get_sum(a,ar);
    while(q--){
        ll l,r;
        fastio::input(l,r);
        fastio::print(ar[r+1]-ar[l]);
    }
}
