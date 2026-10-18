#pragma once
#include "../../template/template.hpp"
namespace tee{
ll Euler(ll N) {
    ll res = N;
    for (ll p = 2; p * p <= N; ++p) {
        if (N % p == 0) {
            res -= res / p;
            while (N % p == 0) N /= p;
        }
    }
    if (N != 1) res -= res / N;
    return res;
}
}