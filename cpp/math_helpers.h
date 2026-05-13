#ifndef MATHHELPERS_H
#define MATHHELPERS_H

namespace constants {
    inline constexpr double PI = 3.141592653589793238462643383279502884;
}

namespace math {
    inline double sinc(double x) {
        if (std::abs(x) < 1e-10) {
            return 1. - x*x/3.0;
        }
        else {
        return std::sin(x)/x;
        }
    }
    inline double sinc2(double x) {
        if (std::abs(x) < 1e-10) {
            return x/3.0;
        }
        else {
            return (std::sin(x) - x*std::cos(x)) / (x*x);
        }
    }
}

#endif