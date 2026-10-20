// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/shortest_path
#include "../template/template.hpp"
#include "../graph/dijkstra.hpp"
signed main(){
    ll n,m,s,t;
    fastio::input(n,m,s,t);
    vvpll g(n);
    rep(i,m){
        ll a,b,c;
        fastio::input(a,b,c);
        g[a].emplace_back(b,c);
        g[b].emplace_back(a,c);
    }
    vpll dist;
    dijkstra(g,s,dist);
    if(dist[t].first==INF){
        fastio::print(-1);
        return 0;
    }
    vll res={t};
    while(res.back()!=s){
        ll v=res.back();
        res.emplace_back(dist[v].second);
    }
    reverse(all(res));
    fastio::print(res);
    return 0;
}
