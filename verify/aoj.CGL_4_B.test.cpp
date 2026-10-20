// competitive-verifier: PROBLEM https://onlinejudge.u-aizu.ac.jp/courses/library/4/CGL/4/CGL_4_B
// competitive-verifier: ERROR 0.000001
#include "../template/template.hpp"
#include "../math/geometry/geometry.hpp"
signed main(){
    ll n;
    fastio::input(n);
    vector<pll> p(n);
    rep(i,n)fastio::input(p[i].first,p[i].second);
    cout<<fixed<<setprecision(10)<<sqrtl((ld)diameter2(p))<<"\n";
}
