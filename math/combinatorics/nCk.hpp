#pragma once
#include "../../template/template.hpp"
namespace tee{
vll fact(1),fact_inv(1),inv(1);
void init_nCk(int SIZE,int Mod) {
    fact.resize(SIZE + 5);
    fact_inv.resize(SIZE + 5);
    inv.resize(SIZE + 5);
    fact[0] = fact[1] = 1;
    fact_inv[0] = fact_inv[1] = 1;
    inv[1] = 1;
    for (int i = 2; i < SIZE + 5; i++) {
        fact[i] = fact[i - 1] * i % Mod;
        inv[i] = Mod - inv[Mod % i] * (Mod / i) % Mod;
        fact_inv[i] = fact_inv[i - 1] * inv[i] % Mod;
    }
    inv[0]=Mod;
}
long long nCk(int n, int k) {
    ll Mod=inv[0];
    if(n < k) return 0;
    if(n < 0 || k < 0) return 0;
    return fact[n] * (fact_inv[k] * fact_inv[n - k] % Mod) % Mod;
}
}