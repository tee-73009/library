#pragma once
#include "../../template/template.hpp"
namespace tee{
vll divisor(ll n){
    vll res;
    for(ll i = 1 ; i*i <= n ; ++i){
        if(n%i == 0){
            res.pb(i);
            if(i*i != n){
                res.pb(n/i);
            }
        }
    }
    sort(all(res));
    return res;
}
}