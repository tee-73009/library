#define PROBLEM "https://judge.yosupo.jp/problem/primality_test"
#include "../template/template.hpp"
#include "library/math/number-theory/is_prime.hpp"
signed main() {
    inl(t);
    rep(i,t){
        ll a;
        fastio::input(a);
        YesNo(is_prime(a));
    }
}