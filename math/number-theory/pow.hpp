#pragma once
#include "../../template/template.hpp"
namespace tee{
template <typename T>
T mypow(T a, ll b) {
    T res = 1;
    while(b > 0) {
        if(b & 1) res = res * a;
        a = a * a;
        b >>= 1;
    }
    return res;
}
ll modpow(ll a, ll b, ll Mod) {
    a = a % Mod;
    if(a<0) a += Mod;
    ll res = 1 % Mod;
    while(b > 0) {
        if(b & 1) res = res * a % Mod;
        a = a * a % Mod;
        b >>= 1;
    }
    return res;
}
}