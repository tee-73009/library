#pragma once
#include "../template/template.hpp"
namespace tee {
template <class T = long long>
struct WeightedUnionFind {
    vector<long long> par, rank, siz;
    vector<T> wt;
    WeightedUnionFind(int n) : par(n,-1), rank(n,0), siz(n,1), wt(n,T()) { }
    int root(int x) {
        if (par[x]==-1) return x;
        int p = par[x];
        int r = root(p);
        wt[x] += wt[p];
        return par[x] = r;
    }
    bool same(int x, int y) {
        return root(x)==root(y);
    }
    T weight(int x) {
        root(x);
        return wt[x];
    }
    T diff(int x, int y) {
        assert(same(x, y));
        return weight(y) - weight(x);
    }
    bool unite(int x, int y, T w) {
        int rx = root(x), ry = root(y);
        if (rx==ry) return false;
        T d = w + weight(x) - weight(y);
        if (rank[rx]<rank[ry]) swap(rx, ry), d = -d;
        par[ry] = rx;
        wt[ry] = d;
        if (rank[rx]==rank[ry]) rank[rx]++;
        siz[rx] += siz[ry];
        return true;
    }
    int size(int x) {
        return siz[root(x)];
    }
};
}