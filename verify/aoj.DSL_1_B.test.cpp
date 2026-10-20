#define PROBLEM "https://onlinejudge.u-aizu.ac.jp/courses/library/3/DSL/1/DSL_1_B"
#include "../template/template.hpp"
#include "../data_structure/weighted_unionfind.hpp"
signed main(){
    ll n,q;
    fastio::input(n,q);
    WeightedUnionFind<ll> uf(n);
    while(q--){
        ll t,x,y;
        fastio::input(t,x,y);
        if(t==0){
            ll z;
            fastio::input(z);
            uf.unite(x,y,z);
        }else{
            if(uf.same(x,y))fastio::print(uf.diff(x,y));
            else fastio::print('?');
        }
    }
}
