// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/factorize
#include "../template/template.hpp"
#include "../math/number-theory/prime_factorize.hpp"
signed main() {
    ll t;
    fastio::input(t);
    while(t--){
        ll a;
        fastio::input(a);
        auto ans=prime_factorize(a);
        fastio::print(ans.size(),ans);
    }
}