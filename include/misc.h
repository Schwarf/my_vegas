//
// Created by andreas on 16.10.24.
//

#ifndef ABS_VEGAS_MISC_H
#define ABS_VEGAS_MISC_H
#include <iostream>
#include <limits>

constexpr double pow_constexpr_double(double base, unsigned int exp) {
    double result = 1;
    while (exp > 0) {
        if (exp % 2 == 1) {
            result *= base;
        }
        base *= base;
        exp /= 2;
    }
    return result;
}

constexpr long long unsigned pow_constexpr(long long unsigned base, unsigned int exp) {
    long long unsigned result = 1;
    while (exp > 0) {
        if (exp % 2 == 1) {
            result *= base;
        }
        base *= base;
        exp /= 2;
    }
    return result;
}

constexpr long long unsigned nth_root(int x, int n) {
    if (x == 0) return 0;
    long long unsigned y = x / 2 > 0 ? x / 2 : 1;
    while (true) {
        long long unsigned y_next = ((n - 1) * y + x / pow_constexpr(y, n - 1)) / n;
        if (y_next == y) {
            break;
        }
        y = y_next;
    }
    return y;
}
#endif //ABS_VEGAS_MISC_H
