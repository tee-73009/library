#pragma once
#include "../template/template.hpp"
namespace tee {
void dijkstra(const vector<vector<pll>> &G, int s, vpll &dis) {
    pqr<pll> pq;
    dis.assign(G.size(), {INF,-1});
    dis[s].first = 0;
    pq.emplace(0, s);
    while (!pq.empty()) {
        auto [d,v]=pq.top();pq.pop();
        if (d > dis[v].fi) continue;
        for (auto [e,c] : G[v]) {
            ll dp = dis[v].fi + c;
            if(dis[e].fi>dp) {
                dis[e] = {dp, v};
                pq.emplace(dp, e);

            }
        }
    }
}
}