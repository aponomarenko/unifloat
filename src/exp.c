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

#include <stddef.h>

#include "unifloat/unifloat.h"
#include "unifloat/exp.h"
#include "unifloat/trig.h"
#include "unifloat/debug.h"

//////////////////////////////////////////////////////////////////////////
//                      exp_UF                                          //
//////////////////////////////////////////////////////////////////////////
Unifloat* exp_UF(Unifloat* x)
{
    Unifloat *res;
    res = mul_UF(x, Log2_e);
    res = call1_arg1(&exp2_UF, res);

    return res;
}

//////////////////////////////////////////////////////////////////////////
//                      exp2_UF                                         //
//////////////////////////////////////////////////////////////////////////
Unifloat* exp2_UF(Unifloat* x)
{
    /*
     * var - the same thing as x, used to not change x in calculations
     * inaccuracy - evaluation of inaccuracy of chain fraction 
     * curNumer, preNumer, curDenom, preDenom - current and previous numerators
     * and denominators of chain fraction
     * a, b_even, b_odd - factors of chain fraction
     * temp - auxiliary variable
     */
    uint k;
    uint i;
    int exponent = 0;
    Unifloat *res, *temp, *inaccuracy, *var, *optmp, *optmp2;
    Unifloat *a, *b_even, *b_odd, *preNumer, *curNumer, *preDenom, *curDenom;
    Unifloat *one;

    if((x->kind == Infinity && x->sign == 1) || x->kind == NaN)
        return clone(x);
    
    if(x->kind == Infinity && x->sign ==  -1) 
        return createZero_UF();

    var = clone(x);
    var->sign = 1;

    /* calculation of integer part of exp2 */
    for(i = var->exp, k = 1; ((int)i) > 0; i--)
    {
        exponent += getMant_UF(var, i) * k;
        if(i < MAX_SIZE_UNIFLOAT * sizeof(int) * 8)
            setMant_UF(var, i, 0);

        k <<= 1;
        if(k == 0x1u << (sizeof(int) * 8 - 1))
        {
            if(x->sign == 1)
                return clone(infinity_UF);
            else
                return createZero_UF();
        }

    }
    normalize_UF(var);
    
    one = createOne_UF();
    /* calculation of fractional part of exp2 */
    if(isZero_UF(var) == false)
    {
        var = call2_arg1(&mul_UF, var, Ln_2);

        /* b[] = 1, 1, 2, 3, 2, 5, 2, 7, 2, 9, ... */
        b_odd = createOne_UF();
        b_even = add_UF(one, one);

        /* a[] = -, x, -x, x, -x, x, -x, ... */
        a = clone(var);
        a->sign = 1;

        /* preNumer = b0 */
        preNumer = createOne_UF();
        /* preDenom = 1 */
        preDenom = createOne_UF();
        /* curNumer = b0*b1 + a1 */
        curNumer = add_UF(one, a);
        /* curDenom = b1 */
        curDenom = createOne_UF();

        k = 2;
        a->sign = -a->sign;

        inaccuracy = createOne_UF();
        temp = createZero_UF();
        
        while((inaccuracy->exp - curDenom->exp - preDenom->exp) >= - (int)PRECISION && isZero_UF(inaccuracy) == 0)
        {    
            if(k%2 == 0)
            {
                copy(curNumer, temp);
                optmp = mul_UF(b_even, curNumer);
                optmp2 = mul_UF(a, preNumer);
                delete_UF(curNumer);
                curNumer = add_UF(optmp, optmp2);
                delete_UFs(optmp, optmp2, NULL);
                copy(temp, preNumer);

                copy(curDenom, temp);
                optmp = mul_UF(b_even, curDenom);
                optmp2 = mul_UF(a, preDenom);
                delete_UF(curDenom);
                curDenom = add_UF(optmp, optmp2);
                delete_UFs(optmp, optmp2, NULL);
                copy(temp, preDenom);

                inaccuracy = call2_arg1(&mul_UF, inaccuracy, b_even);
                b_odd = call2_arg1(&add_UF,b_odd, b_even);
                a->sign =  -a->sign;
            }
            else
            {
                copy(curNumer, temp);

                optmp = mul_UF(b_odd, curNumer);
                optmp2 = mul_UF(a, preNumer);
                delete_UF(curNumer);
                curNumer = add_UF(optmp, optmp2);
                delete_UFs(optmp, optmp2, NULL);

                copy(temp, preNumer);

                copy(curDenom, temp);

                optmp = mul_UF(b_odd, curDenom);
                optmp2 = mul_UF(a, preDenom);
                delete_UF(curDenom);
                curDenom = add_UF(optmp, optmp2);
                delete_UFs(optmp, optmp2, NULL);
                copy(temp, preDenom);

                inaccuracy = call2_arg1(&mul_UF, inaccuracy, b_odd);
                a->sign = -a->sign;
            }
            k++;
        }
        res = div_UF(curNumer, curDenom);
        delete_UFs(curDenom, curNumer, inaccuracy, preDenom, preNumer, b_even, b_odd, temp, a, NULL);
    }
    else
        res = createOne_UF();

    res->exp += exponent;

    if(x->sign == -1)
        res = call2_arg2(&div_UF, one, res);

    delete_UFs(one, var, NULL);

    return res;
}

//////////////////////////////////////////////////////////////////////////
//                      log_UF                                          //
//////////////////////////////////////////////////////////////////////////
Unifloat* log_UF(Unifloat* x)
{
    /*
     * var - the same thing as x, used to not change x in calculations
     * inaccuracy - inaccuracy in newton method
     * exponent - Unifloat number equals var->exp
     * temp - auxiliary variable
     */

    Unifloat *res, *exponent, *var, *optmp, *optmp2;
    Unifloat *one;
   
    int sign;

    if((isInfinity_UF(x) && x->sign == 1) || isNan_UF(x))
        return x;

    if(isZero_UF(x))
    {
        optmp = clone(infinity_UF);
        optmp2 = changeSign_UF(optmp, -1);
        delete_UF(optmp);
        optmp = optmp2;
        return optmp;
    }
    if(x->sign == -1)
        return clone(nan_UF);

    one = createOne_UF();
    if(compare_UF(x, one) == -1)
    {
        var = div_UF(one, x);
        sign = -1;
    }
    else
    {
        var = clone(x);
        sign = 1;
    }
    exponent = convertInteger_UF(var->exp - 1);
    exponent = call2_arg1(&mul_UF, exponent, Ln_2);
    var->exp = 1;

    optmp = sub_UF(var, one);
    optmp = call1_arg1(&log1p_UF,optmp);
    res = add_UF(optmp, exponent);
    res->sign = sign;

    delete_UFs(optmp, var, exponent, one, NULL);

    return res;
}

//////////////////////////////////////////////////////////////////////////
//                      log1p_UF                                        //
//////////////////////////////////////////////////////////////////////////
Unifloat* log1p_UF(Unifloat* x)
{
    /*
     * var - the same thing as x, used to not change x in calculations
     * inaccuracy - evaluation of inaccuracy of chain fraction
     *
     * curNumer, preNumer, curDenom, preDenom - current and previous numerators
     * and denominators of chain fraction
     *
     * a, b - factors of chain fraction
     * temp - auxiliary variable
     */

    int i, k, r;
    Unifloat *res, *temp, *inaccuracy;
    Unifloat *var;
    Unifloat *one;
    Unifloat *exponent;
    Unifloat *inacAddition;
    Unifloat *Zero;
    Unifloat *a, *b, *preNumer, *curNumer, *preDenom, *curDenom, *optmp, *optmp2;

    int sign;

    if((isInfinity_UF(x) && x->sign == 1) || isNan_UF(x) || isZero_UF(x))
        return clone(x);

    optmp = createOne_UF();
    optmp2 = changeSign_UF(optmp, -1);
    r = compare_UF(x, optmp2);
    delete_UFs(optmp, optmp2, NULL);
    if(r == 0)
    {
        optmp = clone(infinity_UF);
        optmp2 = changeSign_UF(optmp, -1);
        delete_UF(optmp);
        return optmp2;
    }

    if(r == -1)
        return clone(nan_UF);

    one = createOne_UF();
    inacAddition = createZero_UF();
    Zero = createZero_UF();
    if(compare_UF(x, Zero) == -1)
    {
        var = add_UF(x, one);
    
        delete_UF(inacAddition);
        inacAddition = clone(x);

        for(i = 1; i <= PRECISION + x->exp; i++)
            setMant_UF(inacAddition, i, 0);

        inacAddition = call1_arg1(&normalize_UF, inacAddition);
        inacAddition = call2_arg1(&div_UF, inacAddition, var);

        var = call2_arg2(&div_UF, one, var);
        var = call2_arg1(&sub_UF, var, one);
        sign = -1;
    }
    else
    {
        var = clone(x);
        sign = 1;
    }

    delete_UF(Zero);

    exponent = createZero_UF();
    /* desire: log(1 + var), 0 <= var < 1 */
    if(var->exp > 1)
    {
        var = call2_arg1(&add_UF, var, one);
        delete_UF(exponent);
        exponent = convertInteger_UF(var->exp - 1);
        exponent = call2_arg1(&mul_UF, exponent, Ln_2);
        var->exp = 1;

        var = call2_arg1(&sub_UF, var, one);
    }
    /* b[] = 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, ...*/    
    b = createOne_UF();

    /* a[] = -, x, 1^2*x, 1^2*x, 2^2*x, 2^2*x, 3^2*x, 3^2*x, ...*/
    a = clone(var);
    
    /* preNumer = b0 */
    preNumer = createZero_UF();
    /* preDenom = 1 */
    preDenom = createOne_UF();
    /* curNumer = b0*b1 + a1 */
    curNumer = clone(a);
    /* curDenom = b1 */
    curDenom = createOne_UF();
    
    k = 2;
    b = call2_arg1(&add_UF, b, one);

    inaccuracy = createOne_UF();
    temp = createZero_UF();

    while((inaccuracy->exp - curDenom->exp - preDenom->exp) >= - (int)PRECISION - 10 && isZero_UF(inaccuracy) == 0)
    {
        copy(curNumer, temp);
        optmp = mul_UF(b, curNumer);
        optmp2 = mul_UF(a, preNumer);
        delete_UF(curNumer);
        curNumer = add_UF(optmp, optmp2);
        delete_UFs(optmp, optmp2, NULL);
        copy(temp, preNumer);

        copy(curDenom, temp);
        optmp = mul_UF(b, curDenom);
        optmp2 = mul_UF(a, preDenom);
        delete_UF(curDenom);
        curDenom = add_UF(optmp, optmp2);
        delete_UFs(optmp,optmp2, NULL);
        copy(temp, preDenom);

        inaccuracy = call2_arg1(&mul_UF, inaccuracy, b);
        k++;

        b = call2_arg1(&add_UF,b, one);
        optmp = convertInteger_UF((k / 2) * (k / 2));
        delete_UF(a);
        a = mul_UF(var, optmp);
        delete_UF(optmp);
    }

    res = div_UF(curNumer, curDenom);
    res = call2_arg1(&add_UF, res, exponent);

    res->sign = sign;
    res = call2_arg1(&add_UF, res, inacAddition);

    delete_UFs(a, b, curNumer, preNumer, curDenom, preDenom, inaccuracy, inacAddition, exponent, one, var, temp, NULL);
    return res;
}

//////////////////////////////////////////////////////////////////////////
//                              sqrt_UF                                 //
//////////////////////////////////////////////////////////////////////////
Unifloat* sqrt_UF(Unifloat* x)
{
    Unifloat *a, *inaccuracy, *res, *optmp;
    int exp;

    if((isInfinity_UF(x) && x->sign == 1) || isNan_UF(x) || isZero_UF(x))
        return x;

    if(x->sign == -1)
        return clone(nan_UF);

    a = clone(x);

    if(x->exp % 2 == 1 || x->exp % 2 == -1)
    {
        exp = (x->exp - 1) / 2;
        a->exp = 1;
    }
    else
    {
        exp = (x->exp - 2) / 2;
        a->exp = 2;
    }

    /* initial approximation */
    res = createOne_UF();
    res->exp++;
    inaccuracy = createOne_UF();

    while(inaccuracy->sign == 1 && inaccuracy->exp >= res->exp - (int)PRECISION && isZero_UF(inaccuracy) == 0)
    {
        copy(res, inaccuracy);
        optmp = div_UF(a, res);
        res = call2_arg1(&add_UF, res, optmp);
        delete_UF(optmp);
        res->exp--;
        inaccuracy = call2_arg1(&sub_UF, inaccuracy, res);
    }

    res->exp += exp;
    delete_UFs(inaccuracy, a, NULL);
    return res;
}
