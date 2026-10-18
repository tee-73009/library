#pragma once
#include "../../template/template.hpp"
namespace tee{
template <typename T>
ll modinv(ll a, ll m) {
    if(m == 0) return 0;
    a %= m;
    if(a < 0) a += m;
    ll b = m, u = 1, v = 0;
    while (b) {
        ll t = a / b;
        a -= t * b; swap(a, b);
        u -= t * v; swap(u, v);
    }
    u %= m; 
    if (u < 0) u += m;
    return u;
}
}