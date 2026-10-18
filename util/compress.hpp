#pragma once
#include "../template/template.hpp"
namespace tee {
//座標圧縮 ex.vec={1,6,2,8,6}→vals={1,2,6,8},vec={0,2,1,3,2}
// O(NlogN)
template<class T>
vector<T> compress(vector<T> &vec){
    auto vals = vec;
    sort(vals.begin(), vals.end());
    vals.erase(unique(vals.begin(), vals.end()), vals.end());
    for(int i = 0; i < vec.size(); i++){
        vec[i] = lower_bound(vals.begin(), vals.end(), vec[i]) - vals.begin();
    }
    return vals;
}
}