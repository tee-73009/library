#include "../template/template.hpp"
namespace std{
//座標圧縮
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