#pragma once
#include "is_prime.hpp"
#include "../../template/template.hpp"
#include "../../data_structure/modint_fast.hpp"
#include "../../util/xorshift.hpp"
namespace tee {
long long pollard(long long N) {
    if (N % 2 == 0) return 2;
    if (is_prime(N)) return N;
    
    using mint = fastmodint;
    mint::set_mod(N);
    
    long long step = 0;
    while (true) {
        mint r = randInt();  // random r
        auto f = [&](mint x) -> mint { return x * x + r; };
        mint x = ++step, y = f(x);
        while (true) {
            long long p = gcd((y - x).get(), (uint64_t)N);
            if (p == 0 || p == N) break;
            if (p != 1) return p;
            x = f(x);
            y = f(f(y));
        }
    }
}

vector<long long> prime_factorize(long long N) {
    if (N == 1) return {};
    long long p = pollard(N);
    if (p == N) return {p};
    vector<long long> left = prime_factorize(p);
    vector<long long> right = prime_factorize(N / p);
    left.insert(left.end(), right.begin(), right.end());
    sort(left.begin(), left.end());
    return left;
}

vpll prime_factorize_pair(long long N) {
    vpll res;
    for (long long p : prime_factorize(N)) {
        if (!res.empty() && res.back().first == p) ++res.back().second;
        else res.push_back({p, 1});
    }
    return res;
}
}