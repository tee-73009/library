#define PROBLEM "https://onlinejudge.u-aizu.ac.jp/courses/library/4/CGL/5/CGL_5_A"
#define ERROR "0.000001"
#include "../template/template.hpp"
#include "../math/geometry/geometry.hpp"
signed main(){
    ioinit;
    ll n;
    cin>>n;
    vector<pair<ld,ld>> p(n);
    rep(i,n)cin>>p[i].first>>p[i].second;
    cout<<fixed<<setprecision(10)<<sqrtl(closest_pair2(p))<<"\n";
}
