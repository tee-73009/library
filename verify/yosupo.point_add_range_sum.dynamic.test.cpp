#define PROBLEM "https://judge.yosupo.jp/problem/point_add_range_sum"
#include "../template/template.hpp"
#include "../data_structure/dynamic_segtree.hpp"
ll op(ll a, ll b) { return a + b; }
ll e() { return 0; }
int main() {
    int n, q;
    cin >> n >> q;
    DynamicSegtree<ll, op, e> seg(n);
    for (int i = 0; i < n; i++) {
        ll a;
        cin >> a;
        seg.set(i, a);
    }
    while (q--) {
        int t;
        cin >> t;
        if (t == 0) {
            int p;
            ll x;
            cin >> p >> x;
            seg.set(p, seg.get(p) + x);
        } else {
            int l, r;
            cin >> l >> r;
            cout << seg.prod(l, r) << "\n";
        }
    }
}
