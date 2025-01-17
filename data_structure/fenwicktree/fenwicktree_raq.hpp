#include "../template/template.hpp"
namespace std {
template <class Abel> struct BIT {
    Abel UNITY_SUM = 0;
    vector<Abel> dat[2];
    int N;

    explicit BIT(int n, Abel unity = 0) : UNITY_SUM(unity), N(n) {
        init(n);
    }
    void init(int n) {
        for (int iter = 0; iter < 2; ++iter)
            dat[iter].assign(n + 1, UNITY_SUM);
    }

    inline void sub_add(int p, int a, Abel x) {
        for (int i = a; i < (int)dat[p].size(); i |= i + 1)
            dat[p][i] = dat[p][i] + x;
    }
    inline void add(int a, int b, Abel x) {
        sub_add(0, a, x * (-a));
        sub_add(1, a, x);
        sub_add(0, b, x * b);
        sub_add(1, b, x * (-1));
    }

    inline Abel sub_sum(int p, int a) {
        Abel res = UNITY_SUM;
        for (int i = a - 1; i >= 0; i = (i & (i + 1)) - 1)
            res = res + dat[p][i];
        return res;
    }
    inline Abel sum(int a, int b) {
        return sub_sum(0, b) + sub_sum(1, b) * b - sub_sum(0, a) - sub_sum(1, a) * a;
    }
    
    friend ostream& operator<<(ostream& os, BIT bit) {
        for (int i = 0; i < bit.N; i++) {
            os << bit.sum(i, i + 1) << " ";
        }
        return os;
    }
};
}