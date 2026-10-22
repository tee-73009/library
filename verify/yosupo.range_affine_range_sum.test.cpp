// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/range_affine_range_sum
#include "../template/template.hpp"
#include "../data_structure/lazysegtree.hpp"
using S=pll; 
using F=pll; 
S op(S a,S b){return{(a.first+b.first)%mod,a.second+b.second};}
S e(){return{0,0};}
S mapping(F f,S x){return{(f.first*x.first+f.second*x.second)%mod,x.second};}
F composition(F f,F g){return{f.first*g.first%mod,(f.first*g.second+f.second)%mod};}
F id(){return{1,0};}
signed main(){
    ll n,q;
    fastio::input(n,q);
    vector<S> a(n);
    rep(i,n){
        fastio::input(a[i].first);
        a[i].second=1;
    }
    LazySegtree<S,op,e,F,mapping,composition,id> seg(a);
    while(q--){
        ll t,l,r;
        fastio::input(t,l,r);
        if(t==0){
            ll b,c;
            fastio::input(b,c);
            seg.apply(l,r,F{b,c});
        }else fastio::print(seg.prod(l,r).first);
    }
}
