#define LINE1 #pragma once
#ifdef USE_LONG_DOUBLE
	#define LINE2 #define ABS_FUNC fabsl
	#define LINE3 #define POW_FUNC powl
	#define LINE4 #define SIN_FUNC sinl
	#define LINE5 #define COS_FUNC cosl
	#define LINE6 #define PI_CONST M_PIl
	#define LINE7 typedef long double FloatType;
#else
	#define LINE2 #define ABS_FUNC fabs
	#define LINE3 #define POW_FUNC pow
	#define LINE4 #define SIN_FUNC sin
	#define LINE5 #define COS_FUNC cos
	#define LINE6 #define PI_CONST M_PI
	#define LINE7 typedef double FloatType;
#endif
LINE1
LINE2
LINE3
LINE4
LINE5
LINE6
LINE7
