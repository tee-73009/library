#pragma once
#include "../template/template.hpp"
namespace tee {
void dijkstra(const vector<vector<pll>> &G, int s, vll &dis) {
    pqr<pll> pq;
    dis.assign(G.size(), INF);
    dis[s] = 0;
    pq.emplace(0, s);
    while (!pq.empty()) {
        auto [d,v]=pq.top();pq.pop();
        if (d > dis[v]) continue;
        for (auto [e,c] : G[v]) {
            ll dp = dis[v] + c;
            if(chmin(dis[e],dp)) {
                pq.emplace(dp, e);
            }
        }
    }
}
}