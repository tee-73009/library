#pragma once
#include "../template/template.hpp"
namespace tee {
struct SegtreeBeats {
    using ll = long long;
  public:
    SegtreeBeats() : SegtreeBeats(0) {}
    explicit SegtreeBeats(int n) : SegtreeBeats(vector<ll>(n, 0)) {}
    explicit SegtreeBeats(const vector<ll> &v) : _n(int(v.size())) {
        d.assign(std::max(1, 4 * _n), Node());
        if (_n) build(1, 0, _n, v);
    }
    int size() const { return _n; }
    void set(int p, ll x) {
        assert(0 <= p && p < _n);
        set_(1, 0, _n, p, x);
    }
    ll get(int p) { return sum(p, p + 1); }
    void chmin(int l, int r, ll x) { range(l, r, 0, x); }
    void chmax(int l, int r, ll x) { range(l, r, 1, x); }
    void add(int l, int r, ll x) { range(l, r, 2, x); }
    // [l, r) を x に置き換える
    void update(int l, int r, ll x) {
        range(l, r, 1, x);
        range(l, r, 0, x);
    }
    ll sum(int l, int r) { return query(l, r, 0); }
    ll min(int l, int r) { return query(l, r, 1); }
    ll max(int l, int r) { return query(l, r, 2); }
    vector<ll> to_vector() {
        vector<ll> res(_n);
        for (int i = 0; i < _n; i++) res[i] = get(i);
        return res;
    }

  private:
    static constexpr ll INF = numeric_limits<ll>::max() / 2;
    struct Node {
        ll sum = 0, mx = -INF, mx2 = -INF, mn = INF, mn2 = INF, ladd = 0;
        int mxc = 0, mnc = 0, len = 0;
    };
    int _n;
    vector<Node> d;
    void apply_add(int k, ll x) {
        Node &n = d[k];
        n.sum += x * n.len;
        n.mx += x, n.mn += x;
        if (n.mx2 != -INF) n.mx2 += x;
        if (n.mn2 != INF) n.mn2 += x;
        n.ladd += x;
    }
    // 前提: mx2 < x < mx
    void apply_chmin(int k, ll x) {
        Node &n = d[k];
        n.sum -= (n.mx - x) * n.mxc;
        if (n.mn == n.mx) n.mn = x;
        else if (n.mn2 == n.mx) n.mn2 = x;
        n.mx = x;
    }
    // 前提: mn < x < mn2
    void apply_chmax(int k, ll x) {
        Node &n = d[k];
        n.sum += (x - n.mn) * n.mnc;
        if (n.mx == n.mn) n.mx = x;
        else if (n.mx2 == n.mn) n.mx2 = x;
        n.mn = x;
    }
    void push(int k) {
        Node &n = d[k];
        for (int c : {2 * k, 2 * k + 1}) {
            if (n.ladd) apply_add(c, n.ladd);
            if (d[c].mx > n.mx) apply_chmin(c, n.mx);
            if (d[c].mn < n.mn) apply_chmax(c, n.mn);
        }
        n.ladd = 0;
    }
    void pull(int k) {
        Node &n = d[k];
        const Node &a = d[2 * k], &b = d[2 * k + 1];
        n.sum = a.sum + b.sum;
        n.len = a.len + b.len;
        if (a.mx > b.mx) {
            n.mx = a.mx, n.mxc = a.mxc;
            n.mx2 = std::max(a.mx2, b.mx);
        } else if (a.mx < b.mx) {
            n.mx = b.mx, n.mxc = b.mxc;
            n.mx2 = std::max(a.mx, b.mx2);
        } else {
            n.mx = a.mx, n.mxc = a.mxc + b.mxc;
            n.mx2 = std::max(a.mx2, b.mx2);
        }
        if (a.mn < b.mn) {
            n.mn = a.mn, n.mnc = a.mnc;
            n.mn2 = std::min(a.mn2, b.mn);
        } else if (a.mn > b.mn) {
            n.mn = b.mn, n.mnc = b.mnc;
            n.mn2 = std::min(a.mn, b.mn2);
        } else {
            n.mn = a.mn, n.mnc = a.mnc + b.mnc;
            n.mn2 = std::min(a.mn2, b.mn2);
        }
    }
    void leaf(int k, ll x) {
        Node &n = d[k];
        n.sum = x, n.mx = n.mn = x, n.mx2 = -INF, n.mn2 = INF;
        n.mxc = n.mnc = 1, n.len = 1, n.ladd = 0;
    }
    void build(int k, int a, int b, const vector<ll> &v) {
        if (b - a == 1) return leaf(k, v[a]);
        int m = (a + b) / 2;
        build(2 * k, a, m, v);
        build(2 * k + 1, m, b, v);
        pull(k);
    }
    void set_(int k, int a, int b, int p, ll x) {
        if (b - a == 1) return leaf(k, x);
        push(k);
        int m = (a + b) / 2;
        if (p < m) set_(2 * k, a, m, p, x);
        else set_(2 * k + 1, m, b, p, x);
        pull(k);
    }
    // type: 0 = chmin, 1 = chmax, 2 = add
    void range(int l, int r, int type, ll x) {
        assert(0 <= l && l <= r && r <= _n);
        if (l == r) return;
        range_(1, 0, _n, l, r, type, x);
    }
    void range_(int k, int a, int b, int l, int r, int type, ll x) {
        if (r <= a || b <= l) return;
        Node &n = d[k];
        if (type == 0 && n.mx <= x) return;
        if (type == 1 && n.mn >= x) return;
        if (l <= a && b <= r) {
            if (type == 2) return apply_add(k, x);
            if (type == 0 && n.mx2 < x) return apply_chmin(k, x);
            if (type == 1 && n.mn2 > x) return apply_chmax(k, x);
        }
        push(k);
        int m = (a + b) / 2;
        range_(2 * k, a, m, l, r, type, x);
        range_(2 * k + 1, m, b, l, r, type, x);
        pull(k);
    }
    // type: 0 = sum, 1 = min, 2 = max
    ll query(int l, int r, int type) {
        assert(0 <= l && l <= r && r <= _n);
        ll res = type == 0 ? 0 : (type == 1 ? INF : -INF);
        if (l == r) return res;
        query_(1, 0, _n, l, r, type, res);
        return res;
    }
    void query_(int k, int a, int b, int l, int r, int type, ll &res) {
        if (r <= a || b <= l) return;
        if (l <= a && b <= r) {
            if (type == 0) res += d[k].sum;
            else if (type == 1) res = std::min(res, d[k].mn);
            else res = std::max(res, d[k].mx);
            return;
        }
        push(k);
        int m = (a + b) / 2;
        query_(2 * k, a, m, l, r, type, res);
        query_(2 * k + 1, m, b, l, r, type, res);
    }
};
}