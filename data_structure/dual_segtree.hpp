#pragma once
#include "../template/template.hpp"
namespace tee {
template <class S, class F, S (*mapping)(F, S), F (*composition)(F, F), F (*id)()>
struct DualSegtree {
  public:
    DualSegtree() : DualSegtree(0) {}
    explicit DualSegtree(int n) : DualSegtree(vector<S>(n, S())) {}
    explicit DualSegtree(const vector<S> &v) : _n(int(v.size())), d(v) {
        size = 1, log = 0;
        while (size < _n) size <<= 1, log++;
        lz = vector<F>(2 * size, id());
    }
    // 位置 p の値を x にする
    void set(int p, S x) {
        assert(0 <= p && p < _n);
        p += size;
        for (int i = log; i >= 1; i--) push(p >> i);
        lz[p] = id();
        d[p - size] = x;
    }
    S get(int p) const {
        assert(0 <= p && p < _n);
        p += size;
        F f = lz[p];
        for (int i = 1; i <= log; i++) f = composition(lz[p >> i], f);
        return mapping(f, d[p - size]);
    }
    void apply(int p, F f) {
        assert(0 <= p && p < _n);
        apply(p, p + 1, f);
    }
    void apply(int l, int r, F f) {
        assert(0 <= l && l <= r && r <= _n);
        if (l == r) return;
        l += size, r += size;
        for (int i = log; i >= 1; i--) {
            if (((l >> i) << i) != l) push(l >> i);
            if (((r >> i) << i) != r) push((r - 1) >> i);
        }
        while (l < r) {
            if (l & 1) all_apply_(l++, f);
            if (r & 1) all_apply_(--r, f);
            l >>= 1, r >>= 1;
        }
    }
    void all_apply(F f) { all_apply_(1, f); }
    vector<S> to_vector() const {
        vector<S> res(_n);
        for (int i = 0; i < _n; i++) res[i] = get(i);
        return res;
    }
  private:
    int _n, size, log;
    vector<S> d;
    vector<F> lz;
    void all_apply_(int k, F f) { lz[k] = composition(f, lz[k]); }
    void push(int k) {
        all_apply_(2 * k, lz[k]);
        all_apply_(2 * k + 1, lz[k]);
        lz[k] = id();
    }
};
}
