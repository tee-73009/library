#pragma once
#include "../../template/template.hpp"
namespace tee {
inline long double EPS = 1e-10;
inline long double torad(long double deg) { return deg * PI / 180; }
inline long double todeg(long double ang) { return ang * 180 / PI; }

template <class T>
int sgn(T x) {
    if constexpr (is_integral_v<T>) return (x > 0) - (x < 0);
    else return x > EPS ? 1 : (x < -EPS ? -1 : 0);
}
template <class T>
pair<long double, long double> to_ld(const pair<T, T> &p) { return {(long double)p.first, (long double)p.second}; }

template <class T>
T dot(const pair<T, T> &a, const pair<T, T> &b) { return a.first * b.first + a.second * b.second; }
template <class T>
T cross(const pair<T, T> &a, const pair<T, T> &b) { return a.first * b.second - a.second * b.first; }
template <class T>
T cross(const pair<T, T> &a, const pair<T, T> &b, const pair<T, T> &c) { return cross(b - a, c - a); }
template <class T>
T norm(const pair<T, T> &a) { return dot(a, a); }
template <class T>
long double abs(const pair<T, T> &a) { return sqrtl((long double)norm(a)); }
template <class T>
long double arg(const pair<T, T> &a) {
    if (sgn(a.first) == 0 && sgn(a.second) == 0) return 0;
    return atan2l((long double)a.second, (long double)a.first);
}
template <class T>
bool eq(const pair<T, T> &a, const pair<T, T> &b) {
    if constexpr (is_integral_v<T>) return a == b;
    else return abs(a - b) <= EPS;
}
template <class T>
pair<T, T> rot90(const pair<T, T> &a) { return {-a.second, a.first}; }
template <class T>
pair<long double, long double> rot(const pair<T, T> &a, long double ang) {
    long double c = cosl(ang), s = sinl(ang), x = a.first, y = a.second;
    return {c * x - s * y, s * x + c * y};
}

template <class T>
void arg_sort(vector<pair<T, T>> &v) {
    auto half = [](const pair<T, T> &p) {  // 偏角が (-π, 0] の側: 0, 原点: 1, 偏角が (0, π] の側: 2
        int sx = sgn(p.first), sy = sgn(p.second);
        if (sx == 0 && sy == 0) return 1;
        return (sy < 0 || (sy == 0 && sx > 0)) ? 0 : 2;
    };
    sort(v.begin(), v.end(), [&](const pair<T, T> &p, const pair<T, T> &q) {
        int hp = half(p), hq = half(q);
        if (hp != hq) return hp < hq;
        int c = sgn(cross(p, q));
        return c != 0 ? c > 0 : norm(p) < norm(q);
    });
}

//a, b を結ぶ線に対する c の位置
//COUNTER_CLOCKWISE(1): c が左側、CLOCKWISE(-1): c が右側
//ONLINE_BACK(2): c, a, b の順に一直線、ONLINE_FRONT(-2): a, b, c の順に一直線
//ON_SEGMENT(0): c が線分 ab 上
constexpr int COUNTER_CLOCKWISE = 1, CLOCKWISE = -1, ONLINE_BACK = 2, ONLINE_FRONT = -2, ON_SEGMENT = 0;
template <class T>
int ccw(const pair<T, T> &a, const pair<T, T> &b, const pair<T, T> &c) {
    pair<T, T> u = b - a, v = c - a;
    int cr = sgn(cross(u, v));
    if (cr > 0) return COUNTER_CLOCKWISE;
    if (cr < 0) return CLOCKWISE;
    if (sgn(dot(u, v)) < 0) return ONLINE_BACK;
    if (sgn(norm(u) - norm(v)) < 0) return ONLINE_FRONT;
    return ON_SEGMENT;
}

// 直線 ab ・線分 ab
template <class T>
bool is_parallel(const pair<T, T> &a, const pair<T, T> &b, const pair<T, T> &c, const pair<T, T> &d) {
    return sgn(cross(b - a, d - c)) == 0;
}
template <class T>
bool is_orthogonal(const pair<T, T> &a, const pair<T, T> &b, const pair<T, T> &c, const pair<T, T> &d) {
    return sgn(dot(b - a, d - c)) == 0;
}
template <class T>
pair<long double, long double> projection(const pair<T, T> &a, const pair<T, T> &b, const pair<T, T> &p) {
    auto A = to_ld(a), B = to_ld(b), P = to_ld(p);
    return A + (B - A) * (dot(P - A, B - A) / norm(B - A));
}
template <class T>
pair<long double, long double> reflection(const pair<T, T> &a, const pair<T, T> &b, const pair<T, T> &p) {
    auto pr = projection(a, b, p);
    return pr + (pr - to_ld(p));
}
template <class T>
bool on_segment(const pair<T, T> &a, const pair<T, T> &b, const pair<T, T> &p) {
    return ccw(a, b, p) == ON_SEGMENT;
}
template <class T>
bool intersect(const pair<T, T> &a, const pair<T, T> &b, const pair<T, T> &c, const pair<T, T> &d) {
    if (eq(a, b)) return ccw(c, d, a) == ON_SEGMENT;  // 長さ 0 の線分 (点)
    if (eq(c, d)) return ccw(a, b, c) == ON_SEGMENT;
    return ccw(a, b, c) * ccw(a, b, d) <= 0 && ccw(c, d, a) * ccw(c, d, b) <= 0;
}
// 直線 ab と直線 cd の交点 (平行でないこと)。線分の場合は intersect で確認してから使う
template <class T>
pair<long double, long double> crosspoint(const pair<T, T> &a, const pair<T, T> &b, const pair<T, T> &c, const pair<T, T> &d) {
    auto A = to_ld(a), B = to_ld(b), C = to_ld(c), D = to_ld(d);
    auto u = B - A, v = D - C;
    return A + u * (cross(C - A, v) / cross(u, v));
}
// 点 a, b の距離
template <class T>
long double distance(const pair<T, T> &a, const pair<T, T> &b) { return abs(a - b); }
// 点 p と直線 ab の距離 / 点 p と線分 ab の距離
template <class T>
long double distance_lp(const pair<T, T> &a, const pair<T, T> &b, const pair<T, T> &p) {
    return fabsl((long double)cross(b - a, p - a)) / abs(b - a);
}
template <class T>
long double distance_sp(const pair<T, T> &a, const pair<T, T> &b, const pair<T, T> &p) {
    if (eq(a, b)) return abs(p - a);
    if (sgn(dot(b - a, p - a)) < 0) return abs(p - a);
    if (sgn(dot(a - b, p - b)) < 0) return abs(p - b);
    return distance_lp(a, b, p);
}
// 線分 ab と線分 cd の距離
template <class T>
long double distance_ss(const pair<T, T> &a, const pair<T, T> &b, const pair<T, T> &c, const pair<T, T> &d) {
    if (intersect(a, b, c, d)) return 0;
    return min({distance_sp(a, b, c), distance_sp(a, b, d), distance_sp(c, d, a), distance_sp(c, d, b)});
}

// 多角形 (頂点を順に並べた vector<pair<T, T>>)
template <class T>
T polygon_area2(const vector<pair<T, T>> &p) {
    T res = 0;
    int n = p.size();
    for (int i = 0; i < n; i++) res += cross(p[i], p[(i + 1) % n]);
    return res;
}
template <class T>
long double polygon_area(const vector<pair<T, T>> &p) { return (long double)polygon_area2(p) / 2; }
// 点の包含判定: 0 (外側), 1 (辺上), 2 (内側)
template <class T>
int contains(const vector<pair<T, T>> &poly, const pair<T, T> &p) {
    bool in = false;
    int n = poly.size();
    for (int i = 0; i < n; i++) {
        pair<T, T> a = poly[i] - p, b = poly[(i + 1) % n] - p;
        if (a.second > b.second) swap(a, b);
        if (sgn(a.second) <= 0 && sgn(b.second) > 0 && sgn(cross(a, b)) < 0) in = !in;
        if (sgn(cross(a, b)) == 0 && sgn(dot(a, b)) <= 0) return 1;
    }
    return in ? 2 : 0;
}
// 凸包 反時計回り、x が最小 (同じなら y が最小) の点から O(N log N)
//   keep_collinear = false: 辺上の点を含めない、true: 含める
template <class T>
vector<pair<T, T>> convex_hull(vector<pair<T, T>> ps, bool keep_collinear = false) {
    sort(ps.begin(), ps.end());
    ps.erase(unique(ps.begin(), ps.end(), [](const pair<T, T> &a, const pair<T, T> &b) { return eq(a, b); }), ps.end());
    int n = ps.size();
    if (n <= 2) return ps;
    auto bad = [&](const pair<T, T> &a, const pair<T, T> &b, const pair<T, T> &c) {
        int cr = sgn(cross(b - a, c - a));
        return keep_collinear ? cr < 0 : cr <= 0;
    };
    vector<pair<T, T>> h(2 * n);
    int k = 0;
    for (int i = 0; i < n; i++) {
        while (k >= 2 && bad(h[k - 2], h[k - 1], ps[i])) k--;
        h[k++] = ps[i];
    }
    for (int i = n - 2, t = k + 1; i >= 0; i--) {
        while (k >= t && bad(h[k - 2], h[k - 1], ps[i])) k--;
        h[k++] = ps[i];
    }
    h.resize(k - 1);
    if (keep_collinear) {
        bool all_collinear = true;
        for (int i = 2; i < n; i++) {
            if (sgn(cross(ps[1] - ps[0], ps[i] - ps[0])) != 0) { all_collinear = false; break; }
        }
        if (all_collinear) return ps;  // 全点が一直線のとき
    }
    return h;
}
// 凸多角形 (反時計回り、convex_hull の結果) の直径の 2 乗 (回転カリパス法) O(N)
template <class T>
T convex_diameter2(const vector<pair<T, T>> &h) {
    int n = h.size();
    if (n <= 1) return 0;
    if (n == 2) return norm(h[0] - h[1]);
    int si = 0, sj = 0;
    for (int k = 1; k < n; k++) {
        if (h[k] < h[si]) si = k;
        if (h[sj] < h[k]) sj = k;
    }
    int i = si, j = sj;
    T best = 0;
    do {
        best = max(best, norm(h[i] - h[j]));
        if (sgn(cross(h[(i + 1) % n] - h[i], h[(j + 1) % n] - h[j])) >= 0) j = (j + 1) % n;
        else i = (i + 1) % n;
    } while (i != si || j != sj);
    return best;
}
// 点集合の直径の 2 乗 (最遠点対)
template <class T>
T diameter2(const vector<pair<T, T>> &ps) { return convex_diameter2(convex_hull(ps)); }

// 最近点対の距離の 2 乗 (点が 2 個以上) O(N log N)
template <class T>
T closest_pair2(vector<pair<T, T>> ps) {
    int n = ps.size();
    assert(n >= 2);
    sort(ps.begin(), ps.end());
    vector<pair<T, T>> buf(n), strip;
    auto by_y = [](const pair<T, T> &a, const pair<T, T> &b) { return a.second < b.second; };
    function<T(int, int)> rec = [&](int l, int r) -> T {
        if (r - l <= 1) return numeric_limits<T>::max();
        int m = (l + r) / 2;
        T mx = ps[m].first;
        T d = min(rec(l, m), rec(m, r));
        merge(ps.begin() + l, ps.begin() + m, ps.begin() + m, ps.begin() + r, buf.begin(), by_y);
        copy(buf.begin(), buf.begin() + (r - l), ps.begin() + l);
        strip.clear();
        for (int i = l; i < r; i++) {
            if ((ps[i].first - mx) * (ps[i].first - mx) >= d) continue;
            for (int j = (int)strip.size() - 1; j >= 0; j--) {
                T dy = ps[i].second - strip[j].second;
                if (dy * dy >= d) break;
                d = min(d, norm(ps[i] - strip[j]));
            }
            strip.push_back(ps[i]);
        }
        return d;
    };
    return rec(0, n);
}

// 円 (中心 c と半径 r を引数で渡す)
// 2 円の共通接線の本数 (4: 離れている, 3: 外接, 2: 交わる, 1: 内接, 0: 一方が他方を含む)
template <class T>
int common_tangent_count(const pair<T, T> &c1, long double r1, const pair<T, T> &c2, long double r2) {
    long double d = abs(c1 - c2);
    if (d > r1 + r2 + EPS) return 4;
    if (fabsl(d - (r1 + r2)) <= EPS) return 3;
    if (fabsl(d - fabsl(r1 - r2)) <= EPS) return 1;
    if (d < fabsl(r1 - r2) - EPS) return 0;
    return 2;
}
// 円 (c, r) と直線 ab の交点 (x, y の昇順、接するときは 1 点)
template <class T>
vector<pair<long double, long double>> crosspoint_cl(const pair<T, T> &c, long double r, const pair<T, T> &a, const pair<T, T> &b) {
    auto center = to_ld(c);
    auto pr = projection(a, b, c);
    long double h2 = norm(center - pr);
    if (h2 > r * r + EPS) return {};
    auto e = to_ld(b - a);
    e = e * (1 / abs(e));
    if (fabsl(h2 - r * r) <= EPS) return {pr};
    long double base = sqrtl(r * r - h2);
    vector<pair<long double, long double>> res = {pr - e * base, pr + e * base};
    sort(res.begin(), res.end());
    return res;
}
// 円 (c1, r1) と円 (c2, r2) の交点 (x, y の昇順、接するときは 1 点、交わらない/同一円のときは空)
template <class T>
vector<pair<long double, long double>> crosspoint_cc(const pair<T, T> &c1, long double r1, const pair<T, T> &c2, long double r2) {
    auto ca = to_ld(c1), cb = to_ld(c2);
    long double d = abs(ca - cb);
    if (d <= EPS) return {};
    if (d > r1 + r2 + EPS || d < fabsl(r1 - r2) - EPS) return {};
    long double x = (d * d + r1 * r1 - r2 * r2) / (2 * d);
    long double h2 = r1 * r1 - x * x;
    auto base = (cb - ca) * (1 / d);
    auto p = ca + base * x;
    if (h2 <= EPS) return {p};
    long double h = sqrtl(h2);
    auto perp = rot90(base);
    vector<pair<long double, long double>> res = {p + perp * h, p - perp * h};
    sort(res.begin(), res.end());
    return res;
}
}  // namespace tee
