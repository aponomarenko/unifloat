/*
    Copyright (c) 2006-2007, 2017 Contributors as noted in the AUTHORS file

    This file is part of libunifloat, the Unifloat C library.

    libunifloat is free software; you can redistribute it and/or modify it
    under the terms of the GNU Lesser General Public License (LGPL) as
    published by the Free Software Foundation; either version 2.1 of the
    License, or (at your option) any later version.

    libunifloat is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser
    General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <math.h>

#include "unifloat/debug.h"
#include "unifloat/libunifloat.h"

void test_point(const char* func_name, double (*func)(double), caller_UF func_UF, double x)
{
    double res, res_UF;
    Unifloat* ures, * ures_UF;

    #ifdef DEBUG_ON
    printf("\n");
    PRINT_InUseStart();
    #endif

    printf("%s: ", func_name);

    res = func(x);
    res_UF = call_UF(func_UF, x);
    
    ures = convertDouble_UF(res);
    ures_UF = convertDouble_UF(res_UF);

    if(compare_UF(ures, ures_UF)!=0)
        printf("FAIL\n");
    else
        printf("OK\n");

    printf("X:    %.17g\n", x);
    printf("LIBM: %.17g\n", res);
    printf("UF:   %.17g\n", res_UF);
    printf("\n");

    delete_UFs(ures, ures_UF, NULL);

    #ifdef DEBUG_ON
    printf("CStrings: %d\n", MEM);
    PRINT_InUseEnd();
    printf("\n");
    #endif
}

void test_point_nx(const char* func_name, double (*func)(int, double), caller_UF_nx func_UF, int n, double x)
{
    double res, res_UF;
    Unifloat* ures, * ures_UF;

    #ifdef DEBUG_ON
    printf("\n");
    PRINT_InUseStart();
    #endif

    printf("%s: ", func_name);

    res = func(n, x);
    res_UF = call_UF_nx(func_UF, n, x);
    
    ures = convertDouble_UF(res);
    ures_UF = convertDouble_UF(res_UF);

    if(compare_UF(ures, ures_UF)!=0)
        printf("FAIL\n");
    else
        printf("OK\n");

    printf("N:    %d\n", n);
    printf("X:    %.17g\n", x);
    printf("LIBM: %.17g\n", res);
    printf("UF:   %.17g\n", res_UF);
    printf("\n");

    delete_UFs(ures, ures_UF, NULL);

    #ifdef DEBUG_ON
    printf("CStrings: %d\n", MEM);
    PRINT_InUseEnd();
    printf("\n");
    #endif
}

int main()
{
    initialize_UF();

    test_point("sin",  &sin,  &sin_UF,   1);
    test_point("cos",  &cos,  &cos_UF,   1);
    test_point("asin", &asin, &asin_UF,  0.5);
    test_point("acos", &acos, &acos_UF,  0.5);
    test_point("tan",  &tan,  &tan_UF,   1);
    test_point("atan", &atan, &atan_UF,  1);
    test_point("exp",  &exp,  &exp_UF,   1.5);
    test_point("log",  &log,  &log_UF,   3);

    printf("\n");

    test_point("j0", &j0, &j0_UF, 4);
    test_point("j1", &j1, &j1_UF, 10);
    test_point("y0", &y0, &y0_UF, 10);
    test_point("y1", &y1, &y1_UF, 10);
    
    test_point_nx("jn", &jn, &jn_UF, 5, 5);
    test_point_nx("yn", &yn, &yn_UF, 5, 5);

    printf("\n");

    test_point("gamma",  &gamma,  &gamma_UF,  5);
    test_point("tgamma", &tgamma, &tgamma_UF, 2.5);
    test_point("lgamma", &lgamma, &lgamma_UF, 5);

    finalize_UF();
}
