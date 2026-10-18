#define PROBLEM "https://judge.yosupo.jp/problem/primality_test"
#include "../template/template.hpp"
#include "../math/number-theory/is_prime.hpp"
signed main() {
    ll t;
    fastio::input(t);
    while(t--){
        ll a;
        fastio::input(a);
        fastio::print(is_prime(a)?"Yes":"No");
    }
}