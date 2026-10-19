#pragma once
#include "../template/template.hpp"
namespace tee {
// Persistent Segment Tree (永続セグメント木)
//   一点更新のたびに新しい「版 (version)」ができ、過去の版にも自由にクエリを投げられる
//   ACL の segtree と同じ op, e の形。版は int で表し、初期状態が版 0、
//   set するたびに 1, 2, 3, ... と増える (どの版から set してもよい)
//   メモリは (初期化 2n + 更新回数 * log n) 程度。添字は 0 始まり、区間は半開区間 [l, r)
// 使い方:
//   ll op(ll a, ll b) { return a + b; }   ll e() { return 0; }
//   PersistentSegtree<ll, op, e> seg(vector<ll>{1, 2, 3});   // 版 0
//   int v1 = seg.set(0, 1, 10);       // 版 0 の位置 1 を 10 にした版
//   int v2 = seg.set(v1, 2, 20);      // 版 v1 から更に更新した版
//   seg.prod(v1, 0, 3);  seg.get(v2, 2);  seg.all_prod(0);
template <class S, S (*op)(S, S), S (*e)()>
struct PersistentSegtree {
  public:
    PersistentSegtree() : PersistentSegtree(0) {}
    explicit PersistentSegtree(int n) : PersistentSegtree(vector<S>(n, e())) {}
    explicit PersistentSegtree(const vector<S> &v) : _n(int(v.size())) {
        nodes.emplace_back();  // 0 番は空の節点
        roots.push_back(_n ? build(0, _n, v) : 0);
    }
    int size() const { return _n; }
    // 版の個数
    int version_count() const { return (int)roots.size(); }
    // 版 ver の位置 p を x にした新しい版を作り、その番号を返す
    int set(int ver, int p, S x) {
        assert(0 <= ver && ver < version_count());
        assert(0 <= p && p < _n);
        int r = set_(roots[ver], 0, _n, p, x);
        roots.push_back(r);
        return (int)roots.size() - 1;
    }
    S get(int ver, int p) const {
        assert(0 <= ver && ver < version_count());
        assert(0 <= p && p < _n);
        int t = roots[ver];
        int a = 0, b = _n;
        while (b - a > 1) {
            int m = (a + b) / 2;
            if (p < m) t = nodes[t].l, b = m;
            else t = nodes[t].r, a = m;
        }
        return nodes[t].val;
    }
    S prod(int ver, int l, int r) const {
        assert(0 <= ver && ver < version_count());
        assert(0 <= l && l <= r && r <= _n);
        return prod_(roots[ver], 0, _n, l, r);
    }
    S all_prod(int ver) const {
        assert(0 <= ver && ver < version_count());
        return nodes[roots[ver]].val;
    }
    // 版 ver で、f(prod(l, r)) が true となる最大の r を返す (f(e()) は true であること)
    template <bool (*f)(S)>
    int max_right(int ver, int l) const {
        return max_right(ver, l, [](S x) { return f(x); });
    }
    template <class G>
    int max_right(int ver, int l, G f) const {
        assert(0 <= ver && ver < version_count());
        assert(0 <= l && l <= _n);
        assert(f(e()));
        if (l == _n) return _n;
        S sm = e();
        int res = mr(roots[ver], 0, _n, l, f, sm);
        return res == -1 ? _n : res;
    }
    // 版 ver で、f(prod(l, r)) が true となる最小の l を返す
    template <bool (*f)(S)>
    int min_left(int ver, int r) const {
        return min_left(ver, r, [](S x) { return f(x); });
    }
    template <class G>
    int min_left(int ver, int r, G f) const {
        assert(0 <= ver && ver < version_count());
        assert(0 <= r && r <= _n);
        assert(f(e()));
        if (r == 0) return 0;
        S sm = e();
        int res = ml(roots[ver], 0, _n, r, f, sm);
        return res == -1 ? 0 : res;
    }

  private:
    struct Node {
        int l = 0, r = 0;
        S val = e();
    };
    int _n;
    vector<Node> nodes;
    vector<int> roots;
    int new_node(int l, int r, S val) {
        nodes.push_back(Node{l, r, val});
        return (int)nodes.size() - 1;
    }
    int build(int a, int b, const vector<S> &v) {
        if (b - a == 1) return new_node(0, 0, v[a]);
        int m = (a + b) / 2;
        int l = build(a, m, v);
        int r = build(m, b, v);
        return new_node(l, r, op(nodes[l].val, nodes[r].val));
    }
    int set_(int t, int a, int b, int p, S x) {
        if (b - a == 1) return new_node(0, 0, x);
        int m = (a + b) / 2;
        int l = nodes[t].l, r = nodes[t].r;
        if (p < m) l = set_(l, a, m, p, x);
        else r = set_(r, m, b, p, x);
        return new_node(l, r, op(nodes[l].val, nodes[r].val));
    }
    S prod_(int t, int a, int b, int l, int r) const {
        if (r <= a || b <= l) return e();
        if (l <= a && b <= r) return nodes[t].val;
        int m = (a + b) / 2;
        return op(prod_(nodes[t].l, a, m, l, r), prod_(nodes[t].r, m, b, l, r));
    }
    template <class G>
    int mr(int t, int a, int b, int l, G &f, S &sm) const {
        if (b <= l) return -1;
        if (l <= a) {
            S nxt = op(sm, nodes[t].val);
            if (f(nxt)) { sm = nxt; return -1; }
            if (b - a == 1) return a;
        }
        int m = (a + b) / 2;
        int res = mr(nodes[t].l, a, m, l, f, sm);
        if (res != -1) return res;
        return mr(nodes[t].r, m, b, l, f, sm);
    }
    template <class G>
    int ml(int t, int a, int b, int r, G &f, S &sm) const {
        if (a >= r) return -1;
        if (b <= r) {
            S nxt = op(nodes[t].val, sm);
            if (f(nxt)) { sm = nxt; return -1; }
            if (b - a == 1) return b;
        }
        int m = (a + b) / 2;
        int res = ml(nodes[t].r, m, b, r, f, sm);
        if (res != -1) return res;
        return ml(nodes[t].l, a, m, r, f, sm);
    }
};
}
