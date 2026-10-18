#pragma once
#include "../template/template.hpp"
#include "compress.hpp"
#include "../data_structure/fenwicktree/fenwicktree.hpp"
namespace tee {
template<class T>
long long inversion(vector<T> a){
    long long res = 0;
    vector<T> b = a;
    compress(b);
    int n = b.size();
    BIT<long long> bit(n);
    rep(i, n){
        res += i - bit.sum(b[i] + 1);
        bit.add(b[i], 1);
    }
    return res;
}
}