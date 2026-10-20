// competitive-verifier: PROBLEM https://onlinejudge.u-aizu.ac.jp/courses/library/3/DSL/2/DSL_2_D
#include "../template/template.hpp"
#include "../data_structure/dual_segtree.hpp"
struct F{int x;bool has;};
int mapping(F f,int x){return f.has?f.x:x;}
F composition(F f,F g){return f.has?f:g;}
F id(){return{0,false};}
signed main(){
    ll n,q;
    fastio::input(n,q);
    DualSegtree<int,F,mapping,composition,id> seg(vector<int>(n,INT_MAX));
    while(q--){
        ll t;
        fastio::input(t);
        if(t==0){
            ll s,tt,x;
            fastio::input(s,tt,x);
            seg.apply(s,tt+1,F{(int)x,true});
        }else{
            ll i;
            fastio::input(i);
            fastio::print(seg.get(i));
        }
    }
}
