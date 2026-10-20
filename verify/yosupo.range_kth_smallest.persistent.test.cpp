// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/range_kth_smallest
#include "../template/template.hpp"
#include "../data_structure/persistent_segtree.hpp"
int op(int a,int b){return a+b;}
int e(){return 0;}
signed main(){
    ll n,q;
    fastio::input(n,q);
    vector<int> a(n);
    rep(i,n)fastio::input(a[i]);
    vector<int> s=a;
    sort(s.begin(),s.end());
    s.erase(unique(s.begin(),s.end()),s.end());
    int m=s.size();
    PersistentSegtree<int,op,e> seg(m);  // 版 i = 先頭 i 個を入れた状態
    vector<int> ver={0};
    rep(i,n){
        int k=lower_bound(s.begin(),s.end(),a[i])-s.begin();
        ver.push_back(seg.set(ver.back(),k,seg.get(ver.back(),k)+1));
    }
    while(q--){
        ll l,r,k;
        fastio::input(l,r,k);
        // 版 ver[r] と ver[l] の差で、個数の和が k+1 以上になる値を二分探索
        int lo=0,hi=m;  // 答えは s[lo..hi) の中にある
        while(hi-lo>1){
            int mid=(lo+hi)/2;
            int c=seg.prod(ver[r],lo,mid)-seg.prod(ver[l],lo,mid);
            if(c>k)hi=mid;
            else k-=c,lo=mid;
        }
        fastio::print(s[lo]);
    }
}
