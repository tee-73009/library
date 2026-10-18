#pragma once
#include "../../template/template.hpp"
namespace tee{
vector<bool> prime_table(int N) {
    vector<bool> isprime(N+1, true);
    isprime[0] = false;
    isprime[1] = false;
    for (int p = 2; p <= N; ++p) {
        if (!isprime[p]) continue;
        for (int q = p * 2; q <= N; q += p) {
            isprime[q] = false;
        }
    }
    return isprime;
}
}