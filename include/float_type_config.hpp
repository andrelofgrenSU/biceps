#define HEADER_LINE1 #pragma once
#ifdef USE_LONG_DOUBLE
    #define HEADER_LINE2 #define ABS_FUNC fabsl
    #define HEADER_LINE3 #define SQRT_FUNC sqrtl
    #define HEADER_LINE4 #define POW_FUNC powl
    #define HEADER_LINE5 #define SIN_FUNC sinl
    #define HEADER_LINE6 #define COS_FUNC cosl
    #define HEADER_LINE7 #define PI_CONST M_PIl
    #define HEADER_LINE8 typedef long double FloatType;
#else
    #define HEADER_LINE2 #define ABS_FUNC fabs
    #define HEADER_LINE3 #define SQRT_FUNC sqrt
    #define HEADER_LINE4 #define POW_FUNC pow
    #define HEADER_LINE5 #define SIN_FUNC sin
    #define HEADER_LINE6 #define COS_FUNC cos
    #define HEADER_LINE7 #define PI_CONST M_PI
    #define HEADER_LINE8 typedef double FloatType;
#endif
HEADER_LINE1
HEADER_LINE2
HEADER_LINE3
HEADER_LINE4
HEADER_LINE5
HEADER_LINE6
HEADER_LINE7
HEADER_LINE8
