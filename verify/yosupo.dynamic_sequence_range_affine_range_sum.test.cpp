// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/dynamic_sequence_range_affine_range_sum
#include "../template/template.hpp"
#include "../data_structure/implicit_treap.hpp"
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
    vector<S> a(n);
    rep(i,n){
        fastio::input(a[i].first);
        a[i].second=1;
    }
    ImplicitTreap<S,op,e,F,mapping,composition,id> t(a);
    while(q--){
        ll ty;
        fastio::input(ty);
        if(ty==0){ll i,x;fastio::input(i,x);t.insert(i,{x,1});}
        else if(ty==1){ll i;fastio::input(i);t.erase(i);}
        else if(ty==2){ll l,r;fastio::input(l,r);t.reverse(l,r);}
        else if(ty==3){ll l,r,b,c;fastio::input(l,r,b,c);t.apply(l,r,{b,c});}
        else{ll l,r;fastio::input(l,r);fastio::print(t.prod(l,r).first);}
    }
}
