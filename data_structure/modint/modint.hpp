#include "../template/template.hpp"
namespace std {
template <int Mod>
struct modint {
    long long val;
    constexpr modint(long long v = 0) noexcept : val(v % Mod) {
        if (val < 0) val += Mod;
    }
    constexpr int getmod() const { return Mod; }
    constexpr modint operator-() const noexcept {
        return val ? Mod - val : 0;
    }
    constexpr modint operator+(const modint& r) const noexcept { return modint(*this) += r; }
    constexpr modint operator-(const modint& r) const noexcept { return modint(*this) -= r; }
    constexpr modint operator*(const modint& r) const noexcept { return modint(*this) *= r; }
    constexpr modint operator/(const modint& r) const noexcept { return modint(*this) /= r; }
    constexpr modint& operator+=(const modint& r) noexcept {
        val += r.val;
        if (val >= Mod) val -= Mod;
        return *this;
    }
    constexpr modint& operator-=(const modint& r) noexcept {
        val -= r.val;
        if (val < 0) val += Mod;
        return *this;
    }
    constexpr modint& operator*=(const modint& r) noexcept {
        val = val * r.val % Mod;
        return *this;
    }
    constexpr modint& operator/=(const modint& r) noexcept {
        long long a = r.val, b = Mod, u = 1, v = 0;
        while (b) {
            long long t = a / b;
            a -= t * b, swap(a, b);
            u -= t * v, swap(u, v);
        }
        val = val * u % Mod;
        if (val < 0) val += Mod;
        return *this;
    }
    constexpr bool operator==(const modint& r) const noexcept {
        return this->val == r.val;
    }
    constexpr bool operator!=(const modint& r) const noexcept {
        return this->val != r.val;
    }
    friend constexpr istream& operator>>(istream& is, modint<Mod>& x) noexcept {
        is >> x.val;
        x.val %= Mod;
        if (x.val < 0) x.val += Mod;
        return is;
    }
    friend constexpr ostream& operator<<(ostream& os, const modint<Mod>& x) noexcept {
        return os << x.val;
    }
    friend constexpr modint<Mod> modpow(const modint<Mod>& r, long long n) noexcept {
        if (n == 0) return 1;
        if (n < 0) return modpow(modinv(r), -n);
        auto t = modpow(r, n / 2);
        t = t * t;
        if (n & 1) t = t * r;
        return t;
    }
    friend constexpr modint<Mod> modinv(const modint<Mod>& r) noexcept {
        long long a = r.val, b = Mod, u = 1, v = 0;
        while (b) {
            long long t = a / b;
            a -= t * b, swap(a, b);
            u -= t * v, swap(u, v);
        }
        return modint<Mod>(u);
    }
};
}