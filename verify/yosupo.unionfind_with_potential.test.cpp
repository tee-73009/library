// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/unionfind_with_potential
#include "../template/template.hpp"
#include "../data_structure/weighted_unionfind.hpp"
signed main(){
    ll n,q;
    fastio::input(n,q);
    WeightedUnionFind<mint> uf(n);
    while(q--){
        ll t;
        fastio::input(t);
        if(t==0){
            ll u,v,x;
            fastio::input(u,v,x);
            if(uf.same(u,v)){
                if(uf.diff(u,v)==x){
                    fastio::print(1);
                }
                else fastio::print(0);
            }
            else{
                uf.unite(u,v,x);
                fastio::print(1);
            }
        }
        else{
            ll u,v;
            fastio::input(u,v);
            if(uf.same(u,v)){
                fastio::print(uf.diff(u,v));
            }
            else{
                fastio::print(-1);
            }
        }
    }
}
