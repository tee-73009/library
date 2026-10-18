#pragma once
#include "../template/template.hpp"
namespace tee {
struct RollingHash {
    static const long long base1 = 1007, base2 = 2009;
    static const long long mod1 = 1000000007, mod2 = 1000000009;
    vector<long long> hash1, hash2, power1, power2;
    RollingHash(const string &S) {
        long long n = (long long)S.size();
        hash1.assign(n+1, 0), hash2.assign(n+1, 0);
        power1.assign(n+1, 1), power2.assign(n+1, 1);
        for (long long i = 0; i < n; ++i) {
            hash1[i+1] = (hash1[i] * base1 + S[i]) % mod1;
            hash2[i+1] = (hash2[i] * base2 + S[i]) % mod2;
            power1[i+1] = (power1[i] * base1) % mod1;
            power2[i+1] = (power2[i] * base2) % mod2;
        }
    }
    inline long long get(long long l, long long r) const {
        long long res1 = hash1[r] - hash1[l] * power1[r-l] % mod1;
        if (res1 < 0) res1 += mod1;
        long long res2 = hash2[r] - hash2[l] * power2[r-l] % mod2;
        if (res2 < 0) res2 += mod2;
        return res1 * mod2 + res2;
    }
    inline long long get() const {
        return hash1.back() * mod2 + hash2.back();
    }
    inline long long getLCP(long long a, long long b) const {
        long long len = min((long long)hash1.size()-a, (long long)hash1.size()-b);
        long long low = 0, high = len;
        while (high - low > 1) {
            long long mid = (low + high) >> 1;
            if (get(a, a+mid) != get(b, b+mid)) high = mid;
            else low = mid;
        }
        return low;
    }
    inline long long getLCP(const RollingHash &T, long long a, long long b) const {
        long long len = min((long long)hash1.size()-a, (long long)hash1.size()-b);
        long long low = 0, high = len;
        while (high - low > 1) {
            long long mid = (low + high) >> 1;
            if (get(a, a+mid) != T.get(b, b+mid)) high = mid;
            else low = mid;
        }
        return low;
    }
};
}