#pragma once
#include "../../template/template.hpp"
namespace tee{
ll basetodec(vll x,ll b){
    ll res=0;
    for(auto i:x){
        res = res*b + i;
    }
    return res;
}
vll dectobase(ll x,ll b){
    vll res;
    while(x){
        res.eb(x%b);
        x/=b;
    }
    reverse(all(res));
    return res;
}
vll basetobase(vll x,ll a,ll b){
    ll t=basetodec(x,a);
    vll res=dectobase(t,b);
    return res;
}
}
