// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/range_chmin_chmax_add_range_sum
#include "../template/template.hpp"
#include "../data_structure/segtree_beats.hpp"
int main() {
    int n, q;
    cin >> n >> q;
    vector<ll> a(n);
    for (auto &x : a) cin >> x;
    SegtreeBeats seg(a);
    while (q--) {
        int t, l, r;
        cin >> t >> l >> r;
        if (t == 3) {
            cout << seg.sum(l, r) << "\n";
            continue;
        }
        ll b;
        cin >> b;
        if (t == 0) seg.chmin(l, r, b);
        else if (t == 1) seg.chmax(l, r, b);
        else seg.add(l, r, b);
    }
}
