// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/binomial_coefficient_prime_mod
#include "../template/template.hpp"
#include "../math/combinatorics/nCk.hpp"
signed main() {
    ll t,m;
    fastio::input(t,m);
    init_nCk(10000000,m);
    while(t--){
        ll n,k;
        fastio::input(n,k);
        fastio::print(nCk(n,k));
    }
}