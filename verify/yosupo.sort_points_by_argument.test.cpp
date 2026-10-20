#define PROBLEM "https://judge.yosupo.jp/problem/sort_points_by_argument"
#include "../template/template.hpp"
#include "../geometry/geometry.hpp"
signed main(){
    ll n;
    fastio::input(n);
    vector<pll> v(n);
    rep(i,n)fastio::input(v[i].first,v[i].second);
    arg_sort(v);
    for(auto&p:v)fastio::print(p.first,p.second);
}
