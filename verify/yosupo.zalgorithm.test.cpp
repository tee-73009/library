#define PROBLEM "https://judge.yosupo.jp/problem/zalgorithm"
#include "../template/template.hpp"
#include "../string/rolling_hash.hpp"
signed main(){
    string s;
    fastio::input(s);
    RollingHash rh(s);
    int n=s.size();
    rep(i,n){
        fastio::out<ll>(rh.getLCP(0,i));
        fastio::out<char>(i+1<n?' ':'\n');
    }
}
