#pragma once
#include "../template/template.hpp"
namespace tee {
template <class S, S (*op)(S, S), S (*e)()>
struct DynamicSegtree {
  public:
    DynamicSegtree() : DynamicSegtree(0) {}
    explicit DynamicSegtree(long long n) : _n(n) { nodes.emplace_back(); }
    long long size() const { return _n; }
    void set(long long p, S x) {
        assert(0 <= p && p < _n);
        if (!root) root = new_node();
        set_(root, 0, _n, p, x);
    }
    S get(long long p) const {
        assert(0 <= p && p < _n);
        int t = root;
        long long a = 0, b = _n;
        while (t && b - a > 1) {
            long long m = a + (b - a) / 2;
            if (p < m) t = nodes[t].l, b = m;
            else t = nodes[t].r, a = m;
        }
        return nodes[t].val;
    }
    S prod(long long l, long long r) const {
        assert(0 <= l && l <= r && r <= _n);
        return prod_(root, 0, _n, l, r);
    }
    S all_prod() const { return nodes[root].val; }
    template <bool (*f)(S)>
    long long max_right(long long l) const {
        return max_right(l, [](S x) { return f(x); });
    }
    template <class G>
    long long max_right(long long l, G f) const {
        assert(0 <= l && l <= _n);
        assert(f(e()));
        S sm = e();
        long long res = mr(root, 0, _n, l, f, sm);
        return res == -1 ? _n : res;
    }
    template <bool (*f)(S)>
    long long min_left(long long r) const {
        return min_left(r, [](S x) { return f(x); });
    }
    template <class G>
    long long min_left(long long r, G f) const {
        assert(0 <= r && r <= _n);
        assert(f(e()));
        S sm = e();
        long long res = ml(root, 0, _n, r, f, sm);
        return res == -1 ? 0 : res;
    }

  private:
    struct Node {
        int l = 0, r = 0;
        S val = e();
    };
    long long _n;
    int root = 0;
    vector<Node> nodes;
    int new_node() {
        nodes.emplace_back();
        return (int)nodes.size() - 1;
    }
    void set_(int t, long long a, long long b, long long p, S x) {
        if (b - a == 1) {
            nodes[t].val = x;
            return;
        }
        long long m = a + (b - a) / 2;
        if (p < m) {
            if (!nodes[t].l) { int c = new_node(); nodes[t].l = c; }
            set_(nodes[t].l, a, m, p, x);
        } else {
            if (!nodes[t].r) { int c = new_node(); nodes[t].r = c; }
            set_(nodes[t].r, m, b, p, x);
        }
        nodes[t].val = op(nodes[nodes[t].l].val, nodes[nodes[t].r].val);
    }
    S prod_(int t, long long a, long long b, long long l, long long r) const {
        if (!t || r <= a || b <= l) return e();
        if (l <= a && b <= r) return nodes[t].val;
        long long m = a + (b - a) / 2;
        return op(prod_(nodes[t].l, a, m, l, r), prod_(nodes[t].r, m, b, l, r));
    }
    template <class G>
    long long mr(int t, long long a, long long b, long long l, G &f, S &sm) const {
        if (b <= l || !t) return -1;
        if (l <= a) {
            S nxt = op(sm, nodes[t].val);
            if (f(nxt)) { sm = nxt; return -1; }
            if (b - a == 1) return a;
        }
        long long m = a + (b - a) / 2;
        long long res = mr(nodes[t].l, a, m, l, f, sm);
        if (res != -1) return res;
        return mr(nodes[t].r, m, b, l, f, sm);
    }
    template <class G>
    long long ml(int t, long long a, long long b, long long r, G &f, S &sm) const {
        if (a >= r || !t) return -1;
        if (b <= r) {
            S nxt = op(nodes[t].val, sm);
            if (f(nxt)) { sm = nxt; return -1; }
            if (b - a == 1) return b;
        }
        long long m = a + (b - a) / 2;
        long long res = ml(nodes[t].r, m, b, r, f, sm);
        if (res != -1) return res;
        return ml(nodes[t].l, a, m, r, f, sm);
    }
};
}