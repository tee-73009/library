#pragma once
#include "../template/template.hpp"
#include "../util/xorshift.hpp"
namespace tee {
template <class S, S (*op)(S, S), S (*e)(), class F, S (*mapping)(F, S), F (*composition)(F, F), F (*id)()>
struct ImplicitTreap {
  private:
    struct Node {
        int l = 0, r = 0, sz = 0;
        unsigned pri = 0;
        S val, acc;
        F lazy;
        bool has_lazy = false, rev = false;
    };
    static vector<Node> &P() {
        static vector<Node> pool = [] {
            vector<Node> p(1);  // 0 番は空の木
            p[0].val = e();
            p[0].acc = e();
            p[0].lazy = id();
            return p;
        }();
        return pool;
    }
    static vector<int> &FREE() {
        static vector<int> f;
        return f;
    }
    static int new_node(const S &x) {
        int t;
        if (!FREE().empty()) {
            t = FREE().back();
            FREE().pop_back();
        } else {
            t = (int)P().size();
            P().emplace_back();
        }
        Node &n = P()[t];
        n.l = n.r = 0;
        n.sz = 1;
        n.pri = randInt();
        n.val = x;
        n.acc = x;
        n.lazy = id();
        n.has_lazy = false;
        n.rev = false;
        return t;
    }
    static void apply_node(int t, const F &f) {
        if (!t) return;
        Node &n = P()[t];
        n.val = mapping(f, n.val);
        n.acc = mapping(f, n.acc);
        n.lazy = n.has_lazy ? composition(f, n.lazy) : f;
        n.has_lazy = true;
    }
    static void rev_node(int t) {
        if (t) P()[t].rev ^= 1;
    }
    static void push(int t) {
        Node &n = P()[t];
        if (n.rev) {
            swap(n.l, n.r);
            rev_node(n.l);
            rev_node(n.r);
            n.rev = false;
        }
        if (n.has_lazy) {
            apply_node(n.l, n.lazy);
            apply_node(n.r, n.lazy);
            n.lazy = id();
            n.has_lazy = false;
        }
    }
    static void pull(int t) {
        Node &n = P()[t];
        n.sz = P()[n.l].sz + P()[n.r].sz + 1;
        n.acc = op(op(P()[n.l].acc, n.val), P()[n.r].acc);
    }
    // [0, k) と [k, n) に分ける
    static pair<int, int> split_(int t, int k) {
        if (!t) return {0, 0};
        push(t);
        int lsz = P()[P()[t].l].sz;
        if (k <= lsz) {
            auto [a, b] = split_(P()[t].l, k);
            P()[t].l = b;
            pull(t);
            return {a, t};
        } else {
            auto [a, b] = split_(P()[t].r, k - lsz - 1);
            P()[t].r = a;
            pull(t);
            return {t, b};
        }
    }
    static int merge_(int a, int b) {
        if (!a) return b;
        if (!b) return a;
        if (P()[a].pri > P()[b].pri) {
            push(a);
            int c = merge_(P()[a].r, b);
            P()[a].r = c;
            pull(a);
            return a;
        } else {
            push(b);
            int c = merge_(a, P()[b].l);
            P()[b].l = c;
            pull(b);
            return b;
        }
    }
    static void collect(int t, vector<S> &out) {
        if (!t) return;
        push(t);
        collect(P()[t].l, out);
        out.push_back(P()[t].val);
        collect(P()[t].r, out);
    }
    int root = 0;
    explicit ImplicitTreap(int r, int) : root(r) {}

  public:
    ImplicitTreap() = default;
    explicit ImplicitTreap(const vector<S> &v) {
        for (auto &x : v) push_back(x);
    }
    int size() const { return P()[root].sz; }
    bool empty() const { return root == 0; }

    // 位置 i (0 <= i <= size) に x を挿入する
    void insert(int i, const S &x) {
        assert(0 <= i && i <= size());
        auto [a, b] = split_(root, i);
        root = merge_(merge_(a, new_node(x)), b);
    }
    void push_back(const S &x) { root = merge_(root, new_node(x)); }
    void push_front(const S &x) { root = merge_(new_node(x), root); }
    // 位置 i の要素を削除して、その値を返す
    S erase(int i) {
        assert(0 <= i && i < size());
        auto [a, bc] = split_(root, i);
        auto [b, c] = split_(bc, 1);
        S res = P()[b].val;
        FREE().push_back(b);
        root = merge_(a, c);
        return res;
    }
    S pop_back() { return erase(size() - 1); }
    S pop_front() { return erase(0); }
    S get(int i) {
        assert(0 <= i && i < size());
        int t = root;
        while (true) {
            push(t);
            int lsz = P()[P()[t].l].sz;
            if (i < lsz) t = P()[t].l;
            else if (i == lsz) return P()[t].val;
            else i -= lsz + 1, t = P()[t].r;
        }
    }
    void set(int i, const S &x) {
        assert(0 <= i && i < size());
        auto [a, bc] = split_(root, i);
        auto [b, c] = split_(bc, 1);
        P()[b].val = x;
        P()[b].acc = x;
        root = merge_(merge_(a, b), c);
    }
    // [l, r) の積 (空区間なら e())
    S prod(int l, int r) {
        assert(0 <= l && l <= r && r <= size());
        auto [a, bc] = split_(root, l);
        auto [b, c] = split_(bc, r - l);
        S res = P()[b].acc;
        root = merge_(merge_(a, b), c);
        return res;
    }
    S all_prod() const { return P()[root].acc; }
    // [l, r) の各要素に作用 f を適用する
    void apply(int l, int r, const F &f) {
        assert(0 <= l && l <= r && r <= size());
        auto [a, bc] = split_(root, l);
        auto [b, c] = split_(bc, r - l);
        apply_node(b, f);
        root = merge_(merge_(a, b), c);
    }
    void apply(int i, const F &f) { apply(i, i + 1, f); }
    // [l, r) を反転する (op が可換のとき)
    void reverse(int l, int r) {
        assert(0 <= l && l <= r && r <= size());
        auto [a, bc] = split_(root, l);
        auto [b, c] = split_(bc, r - l);
        rev_node(b);
        root = merge_(merge_(a, b), c);
    }
    // [l, r) を、先頭が m になるように左へ回転する (l <= m <= r)
    void rotate(int l, int m, int r) {
        assert(0 <= l && l <= m && m <= r && r <= size());
        auto [a, bcd] = split_(root, l);
        auto [b, cd] = split_(bcd, m - l);
        auto [c, d] = split_(cd, r - m);
        root = merge_(merge_(merge_(a, c), b), d);
    }
    vector<S> to_vector() {
        vector<S> res;
        res.reserve(size());
        collect(root, res);
        return res;
    }
    // 分割: this を [0, k) と [k, n) に分けて返す (this は空になる)
    pair<ImplicitTreap, ImplicitTreap> split(int k) {
        assert(0 <= k && k <= size());
        auto [a, b] = split_(root, k);
        root = 0;
        return {ImplicitTreap(a, 0), ImplicitTreap(b, 0)};
    }
    // 結合: this の後ろに other をつなげる (other は空になる)
    void merge(ImplicitTreap &other) {
        root = merge_(root, other.root);
        other.root = 0;
    }
};
}