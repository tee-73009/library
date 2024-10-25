#include "../template/template.hpp"
namespace std{
//累積和(l以上r以下→s[r+1]-s[l])
template<typename T,typename U>
void get_sum(vector<T> &a,vector<U> &sum){
    sum.resize(a.size()+1);
    rep(i,(ll)a.size()){
    sum[i+1]=sum[i]+a[i];
    }
    return;
}
template<typename T,typename U>
vector<U> get_sum(vector<T> &a){
    vector<U> sum;
    get_sum(a,sum);
    return sum;
}
}