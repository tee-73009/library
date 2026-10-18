#pragma once
#include "../../template/template.hpp"
namespace tee{
ll extgcd(ll a, ll b, ll& x, ll& y) {
    if (b == 0) { x = 1; y = 0; return a; }
    ll d = extgcd(b, a % b, y, x);
    y -= a / b * x;
    return d;
}
}