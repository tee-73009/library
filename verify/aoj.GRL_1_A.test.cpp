// competitive-verifier: PROBLEM https://onlinejudge.u-aizu.ac.jp/courses/library/5/GRL/1/GRL_1_A
#include "../template/template.hpp"
#include "../graph/dijkstra.hpp"
signed main(){
    ll v,e,r;
    fastio::input(v,e,r);
    vector<vector<pll>> G(v);
    rep(i,e){
        ll s,t,d;
        fastio::input(s,t,d);
        G[s].push_back({t,d});
    }
    vll dis;
    dijkstra(G,r,dis);
    rep(i,v){
        if(dis[i]==INF)fastio::print(string("INF"));
        else fastio::print(dis[i]);
    }
}
