//
// Created by andreas on 16.10.24.
//

#ifndef ABS_VEGAS_MISC_H
#define ABS_VEGAS_MISC_H
#include <iostream>
#include <limits>
#include <boost/math/special_functions/pow.hpp>
#include <boost/math/tools/roots.hpp>



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

template<long long unsigned X, int N>
struct ValueMapper;

template<> struct ValueMapper<10000,2> {
    static constexpr int NumberOfStratifications = 100;
    static constexpr long long NumberOfHyperCubes = 10000;
};
template<> struct ValueMapper<10000,3> {
    static constexpr int NumberOfStratifications = 21;
    static constexpr long long NumberOfHyperCubes = 9261;
};
template<> struct ValueMapper<10000,4> {
    static constexpr int NumberOfStratifications = 10;
    static constexpr long long NumberOfHyperCubes = 10000;
};
template<> struct ValueMapper<10000,5> {
    static constexpr int NumberOfStratifications = 6;
    static constexpr long long NumberOfHyperCubes = 7776;
};
template<> struct ValueMapper<10000,6> {
    static constexpr int NumberOfStratifications = 4;
    static constexpr long long NumberOfHyperCubes = 4096;
};
template<> struct ValueMapper<10000,7> {
    static constexpr int NumberOfStratifications = 3;
    static constexpr long long NumberOfHyperCubes = 2187;
};
template<> struct ValueMapper<10000,8> {
    static constexpr int NumberOfStratifications = 3;
    static constexpr long long NumberOfHyperCubes = 6561;
};
template<> struct ValueMapper<10000,9> {
    static constexpr int NumberOfStratifications = 2;
    static constexpr long long NumberOfHyperCubes = 512;
};
template<> struct ValueMapper<10000,10> {
    static constexpr int NumberOfStratifications = 2;
    static constexpr long long NumberOfHyperCubes = 1024;
};
template<> struct ValueMapper<10000,11> {
    static constexpr int NumberOfStratifications = 2;
    static constexpr long long NumberOfHyperCubes = 2048;
};
template<> struct ValueMapper<10000,12> {
    static constexpr int NumberOfStratifications = 2;
    static constexpr long long NumberOfHyperCubes = 4096;
};
template<> struct ValueMapper<10000,13> {
    static constexpr int NumberOfStratifications = 2;
    static constexpr long long NumberOfHyperCubes = 8192;
};
template<> struct ValueMapper<10000,14> {
    static constexpr int NumberOfStratifications = 1;
    static constexpr long long NumberOfHyperCubes = 1;
};
template<> struct ValueMapper<10000,15> {
    static constexpr int NumberOfStratifications = 1;
    static constexpr long long NumberOfHyperCubes = 1;
};
template<> struct ValueMapper<100000,2> {
    static constexpr int NumberOfStratifications = 316;
    static constexpr long long NumberOfHyperCubes = 99856;
};
template<> struct ValueMapper<100000,3> {
    static constexpr int NumberOfStratifications = 46;
    static constexpr long long NumberOfHyperCubes = 97336;
};
template<> struct ValueMapper<100000,4> {
    static constexpr int NumberOfStratifications = 17;
    static constexpr long long NumberOfHyperCubes = 83521;
};
template<> struct ValueMapper<100000,5> {
    static constexpr int NumberOfStratifications = 10;
    static constexpr long long NumberOfHyperCubes = 100000;
};
template<> struct ValueMapper<100000,6> {
    static constexpr int NumberOfStratifications = 6;
    static constexpr long long NumberOfHyperCubes = 46656;
};
template<> struct ValueMapper<100000,7> {
    static constexpr int NumberOfStratifications = 5;
    static constexpr long long NumberOfHyperCubes = 78125;
};
template<> struct ValueMapper<100000,8> {
    static constexpr int NumberOfStratifications = 4;
    static constexpr long long NumberOfHyperCubes = 65536;
};
template<> struct ValueMapper<100000,9> {
    static constexpr int NumberOfStratifications = 3;
    static constexpr long long NumberOfHyperCubes = 19683;
};
template<> struct ValueMapper<100000,10> {
    static constexpr int NumberOfStratifications = 3;
    static constexpr long long NumberOfHyperCubes = 59049;
};
template<> struct ValueMapper<100000,11> {
    static constexpr int NumberOfStratifications = 2;
    static constexpr long long NumberOfHyperCubes = 2048;
};
template<> struct ValueMapper<100000,12> {
    static constexpr int NumberOfStratifications = 2;
    static constexpr long long NumberOfHyperCubes = 4096;
};
template<> struct ValueMapper<100000,13> {
    static constexpr int NumberOfStratifications = 2;
    static constexpr long long NumberOfHyperCubes = 8192;
};
template<> struct ValueMapper<100000,14> {
    static constexpr int NumberOfStratifications = 2;
    static constexpr long long NumberOfHyperCubes = 16384;
};
template<> struct ValueMapper<100000,15> {
    static constexpr int NumberOfStratifications = 2;
    static constexpr long long NumberOfHyperCubes = 32768;
};
template<> struct ValueMapper<1000000,2> {
    static constexpr int NumberOfStratifications = 1000;
    static constexpr long long NumberOfHyperCubes = 1000000;
};
template<> struct ValueMapper<1000000,3> {
    static constexpr int NumberOfStratifications = 99;
    static constexpr long long NumberOfHyperCubes = 970299;
};
template<> struct ValueMapper<1000000,4> {
    static constexpr int NumberOfStratifications = 31;
    static constexpr long long NumberOfHyperCubes = 923521;
};
template<> struct ValueMapper<1000000,5> {
    static constexpr int NumberOfStratifications = 15;
    static constexpr long long NumberOfHyperCubes = 759375;
};
template<> struct ValueMapper<1000000,6> {
    static constexpr int NumberOfStratifications = 9;
    static constexpr long long NumberOfHyperCubes = 531441;
};
template<> struct ValueMapper<1000000,7> {
    static constexpr int NumberOfStratifications = 7;
    static constexpr long long NumberOfHyperCubes = 823543;
};
template<> struct ValueMapper<1000000,8> {
    static constexpr int NumberOfStratifications = 5;
    static constexpr long long NumberOfHyperCubes = 390625;
};
template<> struct ValueMapper<1000000,9> {
    static constexpr int NumberOfStratifications = 4;
    static constexpr long long NumberOfHyperCubes = 262144;
};
template<> struct ValueMapper<1000000,10> {
    static constexpr int NumberOfStratifications = 3;
    static constexpr long long NumberOfHyperCubes = 59049;
};
template<> struct ValueMapper<1000000,11> {
    static constexpr int NumberOfStratifications = 3;
    static constexpr long long NumberOfHyperCubes = 177147;
};
template<> struct ValueMapper<1000000,12> {
    static constexpr int NumberOfStratifications = 3;
    static constexpr long long NumberOfHyperCubes = 531441;
};
template<> struct ValueMapper<1000000,13> {
    static constexpr int NumberOfStratifications = 2;
    static constexpr long long NumberOfHyperCubes = 8192;
};
template<> struct ValueMapper<1000000,14> {
    static constexpr int NumberOfStratifications = 2;
    static constexpr long long NumberOfHyperCubes = 16384;
};
template<> struct ValueMapper<1000000,15> {
    static constexpr int NumberOfStratifications = 2;
    static constexpr long long NumberOfHyperCubes = 32768;
};
template<> struct ValueMapper<10000000,2> {
    static constexpr int NumberOfStratifications = 3162;
    static constexpr long long NumberOfHyperCubes = 9998244;
};
template<> struct ValueMapper<10000000,3> {
    static constexpr int NumberOfStratifications = 215;
    static constexpr long long NumberOfHyperCubes = 9938375;
};
template<> struct ValueMapper<10000000,4> {
    static constexpr int NumberOfStratifications = 56;
    static constexpr long long NumberOfHyperCubes = 9834496;
};
template<> struct ValueMapper<10000000,5> {
    static constexpr int NumberOfStratifications = 25;
    static constexpr long long NumberOfHyperCubes = 9765625;
};
template<> struct ValueMapper<10000000,6> {
    static constexpr int NumberOfStratifications = 14;
    static constexpr long long NumberOfHyperCubes = 7529536;
};
template<> struct ValueMapper<10000000,7> {
    static constexpr int NumberOfStratifications = 9;
    static constexpr long long NumberOfHyperCubes = 4782969;
};
template<> struct ValueMapper<10000000,8> {
    static constexpr int NumberOfStratifications = 7;
    static constexpr long long NumberOfHyperCubes = 5764801;
};
template<> struct ValueMapper<10000000,9> {
    static constexpr int NumberOfStratifications = 5;
    static constexpr long long NumberOfHyperCubes = 1953125;
};
template<> struct ValueMapper<10000000,10> {
    static constexpr int NumberOfStratifications = 5;
    static constexpr long long NumberOfHyperCubes = 9765625;
};
template<> struct ValueMapper<10000000,11> {
    static constexpr int NumberOfStratifications = 4;
    static constexpr long long NumberOfHyperCubes = 4194304;
};
template<> struct ValueMapper<10000000,12> {
    static constexpr int NumberOfStratifications = 3;
    static constexpr long long NumberOfHyperCubes = 531441;
};
template<> struct ValueMapper<10000000,13> {
    static constexpr int NumberOfStratifications = 3;
    static constexpr long long NumberOfHyperCubes = 1594323;
};
template<> struct ValueMapper<10000000,14> {
    static constexpr int NumberOfStratifications = 3;
    static constexpr long long NumberOfHyperCubes = 4782969;
};
template<> struct ValueMapper<10000000,15> {
    static constexpr int NumberOfStratifications = 2;
    static constexpr long long NumberOfHyperCubes = 32768;
};
template<> struct ValueMapper<100000000,2> {
    static constexpr int NumberOfStratifications = 10000;
    static constexpr long long NumberOfHyperCubes = 100000000;
};
template<> struct ValueMapper<100000000,3> {
    static constexpr int NumberOfStratifications = 464;
    static constexpr long long NumberOfHyperCubes = 99897344;
};
template<> struct ValueMapper<100000000,4> {
    static constexpr int NumberOfStratifications = 100;
    static constexpr long long NumberOfHyperCubes = 100000000;
};
template<> struct ValueMapper<100000000,5> {
    static constexpr int NumberOfStratifications = 39;
    static constexpr long long NumberOfHyperCubes = 90224199;
};
template<> struct ValueMapper<100000000,6> {
    static constexpr int NumberOfStratifications = 21;
    static constexpr long long NumberOfHyperCubes = 85766121;
};
template<> struct ValueMapper<100000000,7> {
    static constexpr int NumberOfStratifications = 13;
    static constexpr long long NumberOfHyperCubes = 62748517;
};
template<> struct ValueMapper<100000000,8> {
    static constexpr int NumberOfStratifications = 10;
    static constexpr long long NumberOfHyperCubes = 100000000;
};
template<> struct ValueMapper<100000000,9> {
    static constexpr int NumberOfStratifications = 7;
    static constexpr long long NumberOfHyperCubes = 40353607;
};
template<> struct ValueMapper<100000000,10> {
    static constexpr int NumberOfStratifications = 6;
    static constexpr long long NumberOfHyperCubes = 60466176;
};
template<> struct ValueMapper<100000000,11> {
    static constexpr int NumberOfStratifications = 5;
    static constexpr long long NumberOfHyperCubes = 48828125;
};
template<> struct ValueMapper<100000000,12> {
    static constexpr int NumberOfStratifications = 4;
    static constexpr long long NumberOfHyperCubes = 16777216;
};
template<> struct ValueMapper<100000000,13> {
    static constexpr int NumberOfStratifications = 4;
    static constexpr long long NumberOfHyperCubes = 67108864;
};
template<> struct ValueMapper<100000000,14> {
    static constexpr int NumberOfStratifications = 3;
    static constexpr long long NumberOfHyperCubes = 4782969;
};
template<> struct ValueMapper<100000000,15> {
    static constexpr int NumberOfStratifications = 3;
    static constexpr long long NumberOfHyperCubes = 14348907;
};
template<> struct ValueMapper<1000000000,2> {
    static constexpr int NumberOfStratifications = 31622;
    static constexpr long long NumberOfHyperCubes = 999950884;
};
template<> struct ValueMapper<1000000000,3> {
    static constexpr int NumberOfStratifications = 999;
    static constexpr long long NumberOfHyperCubes = 997002999;
};
template<> struct ValueMapper<1000000000,4> {
    static constexpr int NumberOfStratifications = 177;
    static constexpr long long NumberOfHyperCubes = 981506241;
};
template<> struct ValueMapper<1000000000,5> {
    static constexpr int NumberOfStratifications = 63;
    static constexpr long long NumberOfHyperCubes = 992436543;
};
template<> struct ValueMapper<1000000000,6> {
    static constexpr int NumberOfStratifications = 31;
    static constexpr long long NumberOfHyperCubes = 887503681;
};
template<> struct ValueMapper<1000000000,7> {
    static constexpr int NumberOfStratifications = 19;
    static constexpr long long NumberOfHyperCubes = 893871739;
};
template<> struct ValueMapper<1000000000,8> {
    static constexpr int NumberOfStratifications = 13;
    static constexpr long long NumberOfHyperCubes = 815730721;
};
template<> struct ValueMapper<1000000000,9> {
    static constexpr int NumberOfStratifications = 9;
    static constexpr long long NumberOfHyperCubes = 387420489;
};
template<> struct ValueMapper<1000000000,10> {
    static constexpr int NumberOfStratifications = 7;
    static constexpr long long NumberOfHyperCubes = 282475249;
};
template<> struct ValueMapper<1000000000,11> {
    static constexpr int NumberOfStratifications = 6;
    static constexpr long long NumberOfHyperCubes = 362797056;
};
template<> struct ValueMapper<1000000000,12> {
    static constexpr int NumberOfStratifications = 5;
    static constexpr long long NumberOfHyperCubes = 244140625;
};
template<> struct ValueMapper<1000000000,13> {
    static constexpr int NumberOfStratifications = 4;
    static constexpr long long NumberOfHyperCubes = 67108864;
};
template<> struct ValueMapper<1000000000,14> {
    static constexpr int NumberOfStratifications = 4;
    static constexpr long long NumberOfHyperCubes = 268435456;
};
template<> struct ValueMapper<1000000000,15> {
    static constexpr int NumberOfStratifications = 3;
    static constexpr long long NumberOfHyperCubes = 14348907;
};
template<> struct ValueMapper<10000000000,2> {
    static constexpr int NumberOfStratifications = 100000;
    static constexpr long long NumberOfHyperCubes = 10000000000;
};
template<> struct ValueMapper<10000000000,3> {
    static constexpr int NumberOfStratifications = 2154;
    static constexpr long long NumberOfHyperCubes = 9993948264;
};
template<> struct ValueMapper<10000000000,4> {
    static constexpr int NumberOfStratifications = 316;
    static constexpr long long NumberOfHyperCubes = 9971220736;
};
template<> struct ValueMapper<10000000000,5> {
    static constexpr int NumberOfStratifications = 100;
    static constexpr long long NumberOfHyperCubes = 10000000000;
};
template<> struct ValueMapper<10000000000,6> {
    static constexpr int NumberOfStratifications = 46;
    static constexpr long long NumberOfHyperCubes = 9474296896;
};
template<> struct ValueMapper<10000000000,7> {
    static constexpr int NumberOfStratifications = 26;
    static constexpr long long NumberOfHyperCubes = 8031810176;
};
template<> struct ValueMapper<10000000000,8> {
    static constexpr int NumberOfStratifications = 17;
    static constexpr long long NumberOfHyperCubes = 6975757441;
};
template<> struct ValueMapper<10000000000,9> {
    static constexpr int NumberOfStratifications = 12;
    static constexpr long long NumberOfHyperCubes = 5159780352;
};
template<> struct ValueMapper<10000000000,10> {
    static constexpr int NumberOfStratifications = 10;
    static constexpr long long NumberOfHyperCubes = 10000000000;
};
template<> struct ValueMapper<10000000000,11> {
    static constexpr int NumberOfStratifications = 8;
    static constexpr long long NumberOfHyperCubes = 8589934592;
};
template<> struct ValueMapper<10000000000,12> {
    static constexpr int NumberOfStratifications = 6;
    static constexpr long long NumberOfHyperCubes = 2176782336;
};
template<> struct ValueMapper<10000000000,13> {
    static constexpr int NumberOfStratifications = 5;
    static constexpr long long NumberOfHyperCubes = 1220703125;
};
template<> struct ValueMapper<10000000000,14> {
    static constexpr int NumberOfStratifications = 5;
    static constexpr long long NumberOfHyperCubes = 6103515625;
};
template<> struct ValueMapper<10000000000,15> {
    static constexpr int NumberOfStratifications = 4;
    static constexpr long long NumberOfHyperCubes = 1073741824;
};
template<> struct ValueMapper<100000000000,2> {
    static constexpr int NumberOfStratifications = 316227;
    static constexpr long long NumberOfHyperCubes = 99999515529;
};
template<> struct ValueMapper<100000000000,3> {
    static constexpr int NumberOfStratifications = 4641;
    static constexpr long long NumberOfHyperCubes = 99961946721;
};
template<> struct ValueMapper<100000000000,4> {
    static constexpr int NumberOfStratifications = 562;
    static constexpr long long NumberOfHyperCubes = 99757432336;
};
template<> struct ValueMapper<100000000000,5> {
    static constexpr int NumberOfStratifications = 158;
    static constexpr long long NumberOfHyperCubes = 98465804768;
};
template<> struct ValueMapper<100000000000,6> {
    static constexpr int NumberOfStratifications = 68;
    static constexpr long long NumberOfHyperCubes = 98867482624;
};
template<> struct ValueMapper<100000000000,7> {
    static constexpr int NumberOfStratifications = 37;
    static constexpr long long NumberOfHyperCubes = 94931877133;
};
template<> struct ValueMapper<100000000000,8> {
    static constexpr int NumberOfStratifications = 23;
    static constexpr long long NumberOfHyperCubes = 78310985281;
};
template<> struct ValueMapper<100000000000,9> {
    static constexpr int NumberOfStratifications = 16;
    static constexpr long long NumberOfHyperCubes = 68719476736;
};
template<> struct ValueMapper<100000000000,10> {
    static constexpr int NumberOfStratifications = 12;
    static constexpr long long NumberOfHyperCubes = 61917364224;
};
template<> struct ValueMapper<100000000000,11> {
    static constexpr int NumberOfStratifications = 10;
    static constexpr long long NumberOfHyperCubes = 100000000000;
};
template<> struct ValueMapper<100000000000,12> {
    static constexpr int NumberOfStratifications = 8;
    static constexpr long long NumberOfHyperCubes = 68719476736;
};
template<> struct ValueMapper<100000000000,13> {
    static constexpr int NumberOfStratifications = 7;
    static constexpr long long NumberOfHyperCubes = 96889010407;
};
template<> struct ValueMapper<100000000000,14> {
    static constexpr int NumberOfStratifications = 6;
    static constexpr long long NumberOfHyperCubes = 78364164096;
};
template<> struct ValueMapper<100000000000,15> {
    static constexpr int NumberOfStratifications = 5;
    static constexpr long long NumberOfHyperCubes = 30517578125;
};
template<> struct ValueMapper<1000000000000,2> {
    static constexpr int NumberOfStratifications = 1000000;
    static constexpr long long NumberOfHyperCubes = 1000000000000;
};
template<> struct ValueMapper<1000000000000,3> {
    static constexpr int NumberOfStratifications = 9999;
    static constexpr long long NumberOfHyperCubes = 999700029999;
};
template<> struct ValueMapper<1000000000000,4> {
    static constexpr int NumberOfStratifications = 1000;
    static constexpr long long NumberOfHyperCubes = 1000000000000;
};
template<> struct ValueMapper<1000000000000,5> {
    static constexpr int NumberOfStratifications = 251;
    static constexpr long long NumberOfHyperCubes = 996250626251;
};
template<> struct ValueMapper<1000000000000,6> {
    static constexpr int NumberOfStratifications = 99;
    static constexpr long long NumberOfHyperCubes = 941480149401;
};
template<> struct ValueMapper<1000000000000,7> {
    static constexpr int NumberOfStratifications = 51;
    static constexpr long long NumberOfHyperCubes = 897410677851;
};
template<> struct ValueMapper<1000000000000,8> {
    static constexpr int NumberOfStratifications = 31;
    static constexpr long long NumberOfHyperCubes = 852891037441;
};
template<> struct ValueMapper<1000000000000,9> {
    static constexpr int NumberOfStratifications = 21;
    static constexpr long long NumberOfHyperCubes = 794280046581;
};
template<> struct ValueMapper<1000000000000,10> {
    static constexpr int NumberOfStratifications = 15;
    static constexpr long long NumberOfHyperCubes = 576650390625;
};
template<> struct ValueMapper<1000000000000,11> {
    static constexpr int NumberOfStratifications = 12;
    static constexpr long long NumberOfHyperCubes = 743008370688;
};
template<> struct ValueMapper<1000000000000,12> {
    static constexpr int NumberOfStratifications = 9;
    static constexpr long long NumberOfHyperCubes = 282429536481;
};
template<> struct ValueMapper<1000000000000,13> {
    static constexpr int NumberOfStratifications = 8;
    static constexpr long long NumberOfHyperCubes = 549755813888;
};
template<> struct ValueMapper<1000000000000,14> {
    static constexpr int NumberOfStratifications = 7;
    static constexpr long long NumberOfHyperCubes = 678223072849;
};
template<> struct ValueMapper<1000000000000,15> {
    static constexpr int NumberOfStratifications = 6;
    static constexpr long long NumberOfHyperCubes = 470184984576;
};
template<> struct ValueMapper<10000000000000,2> {
    static constexpr int NumberOfStratifications = 3162277;
    static constexpr long long NumberOfHyperCubes = 9999995824729;
};
template<> struct ValueMapper<10000000000000,3> {
    static constexpr int NumberOfStratifications = 21544;
    static constexpr long long NumberOfHyperCubes = 9999516957184;
};
template<> struct ValueMapper<10000000000000,4> {
    static constexpr int NumberOfStratifications = 1778;
    static constexpr long long NumberOfHyperCubes = 9993716528656;
};
template<> struct ValueMapper<10000000000000,5> {
    static constexpr int NumberOfStratifications = 398;
    static constexpr long long NumberOfHyperCubes = 9986547231968;
};
template<> struct ValueMapper<10000000000000,6> {
    static constexpr int NumberOfStratifications = 146;
    static constexpr long long NumberOfHyperCubes = 9685390482496;
};
template<> struct ValueMapper<10000000000000,7> {
    static constexpr int NumberOfStratifications = 71;
    static constexpr long long NumberOfHyperCubes = 9095120158391;
};
template<> struct ValueMapper<10000000000000,8> {
    static constexpr int NumberOfStratifications = 42;
    static constexpr long long NumberOfHyperCubes = 9682651996416;
};
template<> struct ValueMapper<10000000000000,9> {
    static constexpr int NumberOfStratifications = 27;
    static constexpr long long NumberOfHyperCubes = 7625597484987;
};
template<> struct ValueMapper<10000000000000,10> {
    static constexpr int NumberOfStratifications = 19;
    static constexpr long long NumberOfHyperCubes = 6131066257801;
};
template<> struct ValueMapper<10000000000000,11> {
    static constexpr int NumberOfStratifications = 15;
    static constexpr long long NumberOfHyperCubes = 8649755859375;
};
template<> struct ValueMapper<10000000000000,12> {
    static constexpr int NumberOfStratifications = 12;
    static constexpr long long NumberOfHyperCubes = 8916100448256;
};
template<> struct ValueMapper<10000000000000,13> {
    static constexpr int NumberOfStratifications = 10;
    static constexpr long long NumberOfHyperCubes = 10000000000000;
};
template<> struct ValueMapper<10000000000000,14> {
    static constexpr int NumberOfStratifications = 8;
    static constexpr long long NumberOfHyperCubes = 4398046511104;
};
template<> struct ValueMapper<10000000000000,15> {
    static constexpr int NumberOfStratifications = 7;
    static constexpr long long NumberOfHyperCubes = 4747561509943;
};
template<> struct ValueMapper<100000000000000,2> {
    static constexpr int NumberOfStratifications = 10000000;
    static constexpr long long NumberOfHyperCubes = 100000000000000;
};
template<> struct ValueMapper<100000000000000,3> {
    static constexpr int NumberOfStratifications = 46415;
    static constexpr long long NumberOfHyperCubes = 99994258523375;
};
template<> struct ValueMapper<100000000000000,4> {
    static constexpr int NumberOfStratifications = 3162;
    static constexpr long long NumberOfHyperCubes = 99964883083536;
};
template<> struct ValueMapper<100000000000000,5> {
    static constexpr int NumberOfStratifications = 630;
    static constexpr long long NumberOfHyperCubes = 99243654300000;
};
template<> struct ValueMapper<100000000000000,6> {
    static constexpr int NumberOfStratifications = 215;
    static constexpr long long NumberOfHyperCubes = 98771297640625;
};
template<> struct ValueMapper<100000000000000,7> {
    static constexpr int NumberOfStratifications = 99;
    static constexpr long long NumberOfHyperCubes = 93206534790699;
};
template<> struct ValueMapper<100000000000000,8> {
    static constexpr int NumberOfStratifications = 56;
    static constexpr long long NumberOfHyperCubes = 96717311574016;
};
template<> struct ValueMapper<100000000000000,9> {
    static constexpr int NumberOfStratifications = 35;
    static constexpr long long NumberOfHyperCubes = 78815638671875;
};
template<> struct ValueMapper<100000000000000,10> {
    static constexpr int NumberOfStratifications = 25;
    static constexpr long long NumberOfHyperCubes = 95367431640625;
};
template<> struct ValueMapper<100000000000000,11> {
    static constexpr int NumberOfStratifications = 18;
    static constexpr long long NumberOfHyperCubes = 64268410079232;
};
template<> struct ValueMapper<100000000000000,12> {
    static constexpr int NumberOfStratifications = 14;
    static constexpr long long NumberOfHyperCubes = 56693912375296;
};
template<> struct ValueMapper<100000000000000,13> {
    static constexpr int NumberOfStratifications = 11;
    static constexpr long long NumberOfHyperCubes = 34522712143931;
};
template<> struct ValueMapper<100000000000000,14> {
    static constexpr int NumberOfStratifications = 9;
    static constexpr long long NumberOfHyperCubes = 22876792454961;
};
template<> struct ValueMapper<100000000000000,15> {
    static constexpr int NumberOfStratifications = 8;
    static constexpr long long NumberOfHyperCubes = 35184372088832;
};
template<> struct ValueMapper<1000000000000000,2> {
    static constexpr int NumberOfStratifications = 31622776;
    static constexpr long long NumberOfHyperCubes = 999999961946176;
};
template<> struct ValueMapper<1000000000000000,3> {
    static constexpr int NumberOfStratifications = 99999;
    static constexpr long long NumberOfHyperCubes = 999970000299999;
};
template<> struct ValueMapper<1000000000000000,4> {
    static constexpr int NumberOfStratifications = 5623;
    static constexpr long long NumberOfHyperCubes = 999706081460641;
};
template<> struct ValueMapper<1000000000000000,5> {
    static constexpr int NumberOfStratifications = 1000;
    static constexpr long long NumberOfHyperCubes = 1000000000000000;
};
template<> struct ValueMapper<1000000000000000,6> {
    static constexpr int NumberOfStratifications = 316;
    static constexpr long long NumberOfHyperCubes = 995686217814016;
};
template<> struct ValueMapper<1000000000000000,7> {
    static constexpr int NumberOfStratifications = 138;
    static constexpr long long NumberOfHyperCubes = 953133216331392;
};
template<> struct ValueMapper<1000000000000000,8> {
    static constexpr int NumberOfStratifications = 74;
    static constexpr long long NumberOfHyperCubes = 899194740203776;
};
template<> struct ValueMapper<1000000000000000,9> {
    static constexpr int NumberOfStratifications = 46;
    static constexpr long long NumberOfHyperCubes = 922190162669056;
};
template<> struct ValueMapper<1000000000000000,10> {
    static constexpr int NumberOfStratifications = 31;
    static constexpr long long NumberOfHyperCubes = 819628286980801;
};
template<> struct ValueMapper<1000000000000000,11> {
    static constexpr int NumberOfStratifications = 23;
    static constexpr long long NumberOfHyperCubes = 952809757913927;
};
template<> struct ValueMapper<1000000000000000,12> {
    static constexpr int NumberOfStratifications = 17;
    static constexpr long long NumberOfHyperCubes = 582622237229761;
};
template<> struct ValueMapper<1000000000000000,13> {
    static constexpr int NumberOfStratifications = 14;
    static constexpr long long NumberOfHyperCubes = 793714773254144;
};
template<> struct ValueMapper<1000000000000000,14> {
    static constexpr int NumberOfStratifications = 11;
    static constexpr long long NumberOfHyperCubes = 379749833583241;
};
template<> struct ValueMapper<1000000000000000,15> {
    static constexpr int NumberOfStratifications = 10;
    static constexpr long long NumberOfHyperCubes = 1000000000000000;
};
template<> struct ValueMapper<10000000000000000,2> {
    static constexpr int NumberOfStratifications = 100000000;
    static constexpr long long NumberOfHyperCubes = 10000000000000000;
};
template<> struct ValueMapper<10000000000000000,3> {
    static constexpr int NumberOfStratifications = 215443;
    static constexpr long long NumberOfHyperCubes = 9999934692543308;
};
template<> struct ValueMapper<10000000000000000,4> {
    static constexpr int NumberOfStratifications = 10000;
    static constexpr long long NumberOfHyperCubes = 10000000000000000;
};
template<> struct ValueMapper<10000000000000000,5> {
    static constexpr int NumberOfStratifications = 1584;
    static constexpr long long NumberOfHyperCubes = 9971853425639424;
};
template<> struct ValueMapper<10000000000000000,6> {
    static constexpr int NumberOfStratifications = 464;
    static constexpr long long NumberOfHyperCubes = 9979479338254336;
};
template<> struct ValueMapper<10000000000000000,7> {
    static constexpr int NumberOfStratifications = 193;
    static constexpr long long NumberOfHyperCubes = 9974730326005056;
};
template<> struct ValueMapper<10000000000000000,8> {
    static constexpr int NumberOfStratifications = 100;
    static constexpr long long NumberOfHyperCubes = 10000000000000000;
};
template<> struct ValueMapper<10000000000000000,9> {
    static constexpr int NumberOfStratifications = 59;
    static constexpr long long NumberOfHyperCubes = 8662995818654939;
};
template<> struct ValueMapper<10000000000000000,10> {
    static constexpr int NumberOfStratifications = 39;
    static constexpr long long NumberOfHyperCubes = 8140406085191601;
};
template<> struct ValueMapper<10000000000000000,11> {
    static constexpr int NumberOfStratifications = 28;
    static constexpr long long NumberOfHyperCubes = 8293509467471872;
};
template<> struct ValueMapper<10000000000000000,12> {
    static constexpr int NumberOfStratifications = 21;
    static constexpr long long NumberOfHyperCubes = 7355827511386641;
};
template<> struct ValueMapper<10000000000000000,13> {
    static constexpr int NumberOfStratifications = 17;
    static constexpr long long NumberOfHyperCubes = 9904578032905938;
};
template<> struct ValueMapper<10000000000000000,14> {
    static constexpr int NumberOfStratifications = 13;
    static constexpr long long NumberOfHyperCubes = 3937376385699289;
};
template<> struct ValueMapper<10000000000000000,15> {
    static constexpr int NumberOfStratifications = 11;
    static constexpr long long NumberOfHyperCubes = 4177248169415651;
};
template<> struct ValueMapper<100000000000000000,2> {
    static constexpr int NumberOfStratifications = 316227766;
    static constexpr long long NumberOfHyperCubes = 99999999989350752;
};
template<> struct ValueMapper<100000000000000000,3> {
    static constexpr int NumberOfStratifications = 464158;
    static constexpr long long NumberOfHyperCubes = 99999429057832320;
};
template<> struct ValueMapper<100000000000000000,4> {
    static constexpr int NumberOfStratifications = 17782;
    static constexpr long long NumberOfHyperCubes = 99982138977826576;
};
template<> struct ValueMapper<100000000000000000,5> {
    static constexpr int NumberOfStratifications = 2511;
    static constexpr long long NumberOfHyperCubes = 99823677120673552;
};
template<> struct ValueMapper<100000000000000000,6> {
    static constexpr int NumberOfStratifications = 681;
    static constexpr long long NumberOfHyperCubes = 99743056266780080;
};
template<> struct ValueMapper<100000000000000000,7> {
    static constexpr int NumberOfStratifications = 268;
    static constexpr long long NumberOfHyperCubes = 99298698941612032;
};
template<> struct ValueMapper<100000000000000000,8> {
    static constexpr int NumberOfStratifications = 133;
    static constexpr long long NumberOfHyperCubes = 97906861202319840;
};
template<> struct ValueMapper<100000000000000000,9> {
    static constexpr int NumberOfStratifications = 77;
    static constexpr long long NumberOfHyperCubes = 95151694449171440;
};
template<> struct ValueMapper<100000000000000000,10> {
    static constexpr int NumberOfStratifications = 50;
    static constexpr long long NumberOfHyperCubes = 97656250000000000;
};
template<> struct ValueMapper<100000000000000000,11> {
    static constexpr int NumberOfStratifications = 35;
    static constexpr long long NumberOfHyperCubes = 96549157373046880;
};
template<> struct ValueMapper<100000000000000000,12> {
    static constexpr int NumberOfStratifications = 26;
    static constexpr long long NumberOfHyperCubes = 95428956661682176;
};
template<> struct ValueMapper<100000000000000000,13> {
    static constexpr int NumberOfStratifications = 20;
    static constexpr long long NumberOfHyperCubes = 81920000000000000;
};
template<> struct ValueMapper<100000000000000000,14> {
    static constexpr int NumberOfStratifications = 16;
    static constexpr long long NumberOfHyperCubes = 72057594037927936;
};
template<> struct ValueMapper<100000000000000000,15> {
    static constexpr int NumberOfStratifications = 13;
    static constexpr long long NumberOfHyperCubes = 51185893014090760;
};
template<> struct ValueMapper<1000000000000000000,2> {
    static constexpr int NumberOfStratifications = 1000000000;
    static constexpr long long NumberOfHyperCubes = 1000000000000000000;
};
template<> struct ValueMapper<1000000000000000000,3> {
    static constexpr int NumberOfStratifications = 999999;
    static constexpr long long NumberOfHyperCubes = 999997000002999936;
};
template<> struct ValueMapper<1000000000000000000,4> {
    static constexpr int NumberOfStratifications = 31622;
    static constexpr long long NumberOfHyperCubes = 999901770412381440;
};
template<> struct ValueMapper<1000000000000000000,5> {
    static constexpr int NumberOfStratifications = 3981;
    static constexpr long long NumberOfHyperCubes = 999909945163943936;
};
template<> struct ValueMapper<1000000000000000000,6> {
    static constexpr int NumberOfStratifications = 999;
    static constexpr long long NumberOfHyperCubes = 994014980014994048;
};
template<> struct ValueMapper<1000000000000000000,7> {
    static constexpr int NumberOfStratifications = 372;
    static constexpr long long NumberOfHyperCubes = 985826706403442688;
};
template<> struct ValueMapper<1000000000000000000,8> {
    static constexpr int NumberOfStratifications = 177;
    static constexpr long long NumberOfHyperCubes = 963354501121950080;
};
template<> struct ValueMapper<1000000000000000000,9> {
    static constexpr int NumberOfStratifications = 99;
    static constexpr long long NumberOfHyperCubes = 913517247483640960;
};
template<> struct ValueMapper<1000000000000000000,10> {
    static constexpr int NumberOfStratifications = 63;
    static constexpr long long NumberOfHyperCubes = 984930291881790848;
};
template<> struct ValueMapper<1000000000000000000,11> {
    static constexpr int NumberOfStratifications = 43;
    static constexpr long long NumberOfHyperCubes = 929293739471222656;
};
template<> struct ValueMapper<1000000000000000000,12> {
    static constexpr int NumberOfStratifications = 31;
    static constexpr long long NumberOfHyperCubes = 787662783788549760;
};
template<> struct ValueMapper<1000000000000000000,13> {
    static constexpr int NumberOfStratifications = 24;
    static constexpr long long NumberOfHyperCubes = 876488338465357824;
};
template<> struct ValueMapper<1000000000000000000,14> {
    static constexpr int NumberOfStratifications = 19;
    static constexpr long long NumberOfHyperCubes = 799006685782884096;
};
template<> struct ValueMapper<1000000000000000000,15> {
    static constexpr int NumberOfStratifications = 15;
    static constexpr long long NumberOfHyperCubes = 437893890380859392;
};


#endif //ABS_VEGAS_MISC_H
