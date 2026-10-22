#pragma once
#include "../template/template.hpp"
namespace tee {
template <class S, S (*op)(S, S), S (*e)()>
struct Segtree {
  public:
    Segtree() : Segtree(0) {}
    explicit Segtree(int n) : Segtree(vector<S>(n, e())) {}
    explicit Segtree(const vector<S> &v) : _n(int(v.size())) {
        sz = 1, log = 0;
        while (sz < _n) sz <<= 1, log++;
        d = vector<S>(2 * sz, e());
        for (int i = 0; i < _n; i++) d[sz + i] = v[i];
        for (int i = sz - 1; i >= 1; i--) update(i);
    }
    int size() const { return _n; }
    void set(int p, S x) {
        assert(0 <= p && p < _n);
        p += sz;
        d[p] = x;
        for (int i = 1; i <= log; i++) update(p >> i);
    }
    S get(int p) const {
        assert(0 <= p && p < _n);
        return d[p + sz];
    }
    S prod(int l, int r) const {
        assert(0 <= l && l <= r && r <= _n);
        S sml = e(), smr = e();
        l += sz, r += sz;
        while (l < r) {
            if (l & 1) sml = op(sml, d[l++]);
            if (r & 1) smr = op(d[--r], smr);
            l >>= 1, r >>= 1;
        }
        return op(sml, smr);
    }
    S all_prod() const { return d[1]; }
    template <bool (*f)(S)>
    int max_right(int l) const {
        return max_right(l, [](S x) { return f(x); });
    }
    template <class G>
    int max_right(int l, G f) const {
        assert(0 <= l && l <= _n);
        assert(f(e()));
        if (l == _n) return _n;
        l += sz;
        S sm = e();
        do {
            while (l % 2 == 0) l >>= 1;
            if (!f(op(sm, d[l]))) {
                while (l < sz) {
                    l = 2 * l;
                    if (f(op(sm, d[l]))) sm = op(sm, d[l++]);
                }
                return l - sz;
            }
            sm = op(sm, d[l++]);
        } while ((l & -l) != l);
        return _n;
    }
    template <bool (*f)(S)>
    int min_left(int r) const {
        return min_left(r, [](S x) { return f(x); });
    }
    template <class G>
    int min_left(int r, G f) const {
        assert(0 <= r && r <= _n);
        assert(f(e()));
        if (r == 0) return 0;
        r += sz;
        S sm = e();
        do {
            r--;
            while (r > 1 && (r % 2)) r >>= 1;
            if (!f(op(d[r], sm))) {
                while (r < sz) {
                    r = 2 * r + 1;
                    if (f(op(d[r], sm))) sm = op(d[r--], sm);
                }
                return r + 1 - sz;
            }
            sm = op(d[r], sm);
        } while ((r & -r) != r);
        return 0;
    }
    vector<S> to_vector() const { return vector<S>(d.begin() + sz, d.begin() + sz + _n); }

  private:
    int _n, sz, log;
    vector<S> d;
    void update(int k) { d[k] = op(d[2 * k], d[2 * k + 1]); }
};
}
