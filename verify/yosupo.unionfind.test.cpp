#define PROBLEM "https://judge.yosupo.jp/problem/unionfind"
#include "../template/template.hpp"
#include "../data_structure/unionfind.hpp"
signed main() {
    ll n,q;
    fastio::input(n,q);
    UnionFind uf(n);
    while(q--){
        ll t, u, v;
        fastio::input(t,u,v);
        if (t == 0) uf.unite(u, v);
        else fastio::print((uf.same(u, v) ? 1 : 0));
    }
}