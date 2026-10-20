// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/many_aplusb_128bit
#include "../template/template.hpp"
signed main(){
    ll t;
    fastio::input(t);
    while(t--){
        __int128_t a,b;
        fastio::input(a,b);
        fastio::print(a+b);
    }
}
