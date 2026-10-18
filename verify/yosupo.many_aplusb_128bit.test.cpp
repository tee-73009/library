#define PROBLEM "https://judge.yosupo.jp/problem/many_aplusb_128bit"
#include "../template/template.hpp"
signed main() {
    inl(t);
    rep(i,t){
        __int128_t a,b;
        input(a,b);
        print(a+b);
    }
}