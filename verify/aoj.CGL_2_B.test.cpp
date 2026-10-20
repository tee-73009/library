#define PROBLEM "https://onlinejudge.u-aizu.ac.jp/courses/library/4/CGL/2/CGL_2_B"
#include "../template/template.hpp"
#include "../math/geometry/geometry.hpp"
signed main(){
    ll q;
    fastio::input(q);
    while(q--){
        pll p0,p1,p2,p3;
        fastio::input(p0.first,p0.second,p1.first,p1.second,p2.first,p2.second,p3.first,p3.second);
        fastio::print(intersect(p0,p1,p2,p3)?1:0);
    }
}
