#pragma once
#include "../../template/template.hpp"
namespace tee{
vll Euler_table(int N) {
    vll euler(N+1);
    iota(all(euler),0);
    for (int i = 2; i <= N; ++i) {
        if (euler[i]!=i) continue;
        for (int j = i; j <= N; j += i) {
            euler[j] -= euler[j] / i;
        }
    }
    return euler;
}
}