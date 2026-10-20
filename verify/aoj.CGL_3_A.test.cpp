#define PROBLEM "https://onlinejudge.u-aizu.ac.jp/courses/library/4/CGL/3/CGL_3_A"
#include "../template/template.hpp"
#include "../math/geometry/geometry.hpp"
signed main(){
    ll n;
    fastio::input(n);
    vector<pll> p(n);
    rep(i,n)fastio::input(p[i].first,p[i].second);
    ll a2=llabs(polygon_area2(p));
    cout<<a2/2<<"."<<(a2%2?5:0)<<"\n";
}
