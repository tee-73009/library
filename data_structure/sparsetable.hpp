#pragma once
#include "../template/template.hpp"
namespace tee {
template <class S, S (*op)(S, S)>
struct SparseTable {
    SparseTable() = default;
    explicit SparseTable(const vector<S>& v) { init(v); }
    void init(const vector<S>& v) {
        n = (int)v.size();
        dat.clear();
        if (n == 0) return;
        dat.push_back(v);
        for (int j = 1; (1 << j) <= n; ++j) {
            int len = n - (1 << j) + 1;
            dat.emplace_back(len);
            for (int i = 0; i < len; ++i) {
                dat[j][i] = op(dat[j - 1][i], dat[j - 1][i + (1 << (j - 1))]);
            }
        }
    }
    S prod(int l, int r) const {
        assert(0 <= l && l < r && r <= n);
        int j = 31 - __builtin_clz(r - l);
        return op(dat[j][l], dat[j][r - (1 << j)]);
    }
    S get(int i) const {
        assert(0 <= i && i < n);
        return dat[0][i];
    }
    int size() const { return n; }
  private:
    int n = 0;
    vector<vector<S>> dat;
};
}