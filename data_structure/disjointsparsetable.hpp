#pragma once
#include "../template/template.hpp"
namespace tee {
template <class S, S (*op)(S, S)>
struct DisjointSparseTable {
    DisjointSparseTable() = default;
    explicit DisjointSparseTable(const vector<S>& v) { init(v); }

    void init(const vector<S>& v) {
        n = (int)v.size();
        int lg = 0;
        while ((1 << lg) < n) ++lg;
        dat.assign(lg, v);
        a = v;
        for (int k = 0; k < lg; ++k) {
            int half = 1 << k;
            for (int mid = half; mid < n; mid += 2 * half) {
                for (int i = mid - 2; i >= mid - half; --i) {
                    dat[k][i] = op(a[i], dat[k][i + 1]);
                }
                int right_end = min(mid + half, n);
                for (int i = mid + 1; i < right_end; ++i) {
                    dat[k][i] = op(dat[k][i - 1], a[i]);
                }
            }
        }
    }
    S prod(int l, int r) const {
        assert(0 <= l && l < r && r <= n);
        if (r - l == 1) return a[l];
        --r;
        int k = 31 - __builtin_clz(l ^ r);
        return op(dat[k][l], dat[k][r]);
    }
    S get(int i) const {
        assert(0 <= i && i < n);
        return a[i];
    }
    int size() const { return n; }
  private:
    int n = 0;
    vector<S> a;
    vector<vector<S>> dat;
};
}