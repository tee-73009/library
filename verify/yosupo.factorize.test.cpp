#define PROBLEM "https://judge.yosupo.jp/problem/factorize"
#include "../template/template.hpp"
#include "../math/number-theory/prime_factorize.hpp"
signed main() {
    inl(t);
    rep(i,t){
        ll a;
        fastio::input(a);
        auto ans=prime_factorize(a);
        print(ans.size(),ans);
    }
}