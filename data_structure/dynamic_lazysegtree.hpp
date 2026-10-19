#pragma once
#include "../template/template.hpp"
namespace tee {
template <class S, S (*op)(S, S), S (*e)(), class F, S (*mapping)(F, S), F (*composition)(F, F), F (*id)()>
struct DynamicLazySegtree {
  public:
    DynamicLazySegtree() : DynamicLazySegtree(0) {}
    explicit DynamicLazySegtree(long long n, S init = e()) : _n(n), leaf(init) {
        nodes.emplace_back();
        nodes[0].val = n > 0 ? pw(n) : e();
    }
    long long size() const { return _n; }
    void set(long long p, S x) {
        assert(0 <= p && p < _n);
        set_(0, 0, _n, p, x);
    }
    S get(long long p) {
        assert(0 <= p && p < _n);
        return prod(p, p + 1);
    }
    S prod(long long l, long long r) {
        assert(0 <= l && l <= r && r <= _n);
        if (l == r) return e();
        return prod_(0, 0, _n, l, r);
    }
    S all_prod() const { return nodes[0].val; }
    void apply(long long p, F f) {
        assert(0 <= p && p < _n);
        apply(p, p + 1, f);
    }
    void apply(long long l, long long r, F f) {
        assert(0 <= l && l <= r && r <= _n);
        if (l == r) return;
        apply_(0, 0, _n, l, r, f);
    }
    void all_apply(F f) { apply_node(0, f); }
    template <bool (*f)(S)>
    long long max_right(long long l) {
        return max_right(l, [](S x) { return f(x); });
    }
    template <class G>
    long long max_right(long long l, G f) {
        assert(0 <= l && l <= _n);
        assert(f(e()));
        if (l == _n) return _n;
        S sm = e();
        long long res = mr(0, 0, _n, l, f, sm);
        return res == -1 ? _n : res;
    }
    template <bool (*f)(S)>
    long long min_left(long long r) {
        return min_left(r, [](S x) { return f(x); });
    }
    template <class G>
    long long min_left(long long r, G f) {
        assert(0 <= r && r <= _n);
        assert(f(e()));
        if (r == 0) return 0;
        S sm = e();
        long long res = ml(0, 0, _n, r, f, sm);
        return res == -1 ? 0 : res;
    }

  private:
    struct Node {
        int l = 0, r = 0;  // 0 は「まだ無い」(根は子になれないので 0 を使える)
        S val = e();
        F lz = id();
    };
    long long _n;
    S leaf;
    vector<Node> nodes;
    map<long long, S> memo;
    // 初期値 leaf を w 個 op で並べたもの
    S pw(long long w) {
        if (w == 1) return leaf;
        auto it = memo.find(w);
        if (it != memo.end()) return it->second;
        S a = pw((w + 1) / 2), b = pw(w / 2);
        return memo[w] = op(a, b);
    }
    void apply_node(int t, F f) {
        nodes[t].val = mapping(f, nodes[t].val);
        nodes[t].lz = composition(f, nodes[t].lz);
    }
    void push(int t, long long a, long long b) {
        long long m = a + (b - a) / 2;
        if (!nodes[t].l) {
            int c = (int)nodes.size();
            nodes.emplace_back();
            nodes[c].val = pw(m - a);
            nodes[t].l = c;
        }
        if (!nodes[t].r) {
            int c = (int)nodes.size();
            nodes.emplace_back();
            nodes[c].val = pw(b - m);
            nodes[t].r = c;
        }
        F f = nodes[t].lz;
        apply_node(nodes[t].l, f);
        apply_node(nodes[t].r, f);
        nodes[t].lz = id();
    }
    void pull(int t) { nodes[t].val = op(nodes[nodes[t].l].val, nodes[nodes[t].r].val); }
    void set_(int t, long long a, long long b, long long p, S x) {
        if (b - a == 1) {
            nodes[t].val = x;
            nodes[t].lz = id();
            return;
        }
        push(t, a, b);
        long long m = a + (b - a) / 2;
        if (p < m) set_(nodes[t].l, a, m, p, x);
        else set_(nodes[t].r, m, b, p, x);
        pull(t);
    }
    S prod_(int t, long long a, long long b, long long l, long long r) {
        if (r <= a || b <= l) return e();
        if (l <= a && b <= r) return nodes[t].val;
        push(t, a, b);
        long long m = a + (b - a) / 2;
        return op(prod_(nodes[t].l, a, m, l, r), prod_(nodes[t].r, m, b, l, r));
    }
    void apply_(int t, long long a, long long b, long long l, long long r, F f) {
        if (r <= a || b <= l) return;
        if (l <= a && b <= r) return apply_node(t, f);
        push(t, a, b);
        long long m = a + (b - a) / 2;
        apply_(nodes[t].l, a, m, l, r, f);
        apply_(nodes[t].r, m, b, l, r, f);
        pull(t);
    }
    template <class G>
    long long mr(int t, long long a, long long b, long long l, G &f, S &sm) {
        if (b <= l) return -1;
        if (l <= a) {
            S nxt = op(sm, nodes[t].val);
            if (f(nxt)) { sm = nxt; return -1; }
            if (b - a == 1) return a;
        }
        push(t, a, b);
        long long m = a + (b - a) / 2;
        long long res = mr(nodes[t].l, a, m, l, f, sm);
        if (res != -1) return res;
        return mr(nodes[t].r, m, b, l, f, sm);
    }
    template <class G>
    long long ml(int t, long long a, long long b, long long r, G &f, S &sm) {
        if (a >= r) return -1;
        if (b <= r) {
            S nxt = op(nodes[t].val, sm);
            if (f(nxt)) { sm = nxt; return -1; }
            if (b - a == 1) return b;
        }
        push(t, a, b);
        long long m = a + (b - a) / 2;
        long long res = ml(nodes[t].r, m, b, r, f, sm);
        if (res != -1) return res;
        return ml(nodes[t].l, a, m, r, f, sm);
    }
};
}
