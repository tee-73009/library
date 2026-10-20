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
    }
    vpll dist;
    dijkstra(g,s,dist);
    if(dist[t].first==INF){
        fastio::print(-1);
        return 0;
    }
    vll res={t};
    while(res.back()!=s)res.emplace_back(dist[res.back()].second);
    reverse(all(res));
    ll k=res.size()-1;
    fastio::print(dist[t].first,k);
    rep(i,k)fastio::print(res[i],res[i+1]);
}