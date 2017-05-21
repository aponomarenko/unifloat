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

#include <stddef.h> // NULL

#include "unifloat/unifloat.h"
#include "unifloat/gamma.h"
#include "unifloat/exp.h"

Unifloat* gamma_UF(Unifloat* x)
{
    Unifloat* res, * Half, * One, * Two, * int_x, * Ten, * z, * new_x, * optmp;
    Unifloat* huge_val = clone(infinity_UF);
    CString* cstr;
    int i, num, error;
    Two = createOne_UF();
    Two->exp++;
    new_x = clone(x);

    if((x->sign == -1) && isNormal_UF(x))
    {
        int_x = clone(x);
        for(i = int_x->exp + 1; i < PRECISION; i++)
            setMant_UF(int_x, i, 0);

        int_x = call2_arg1(&div_UF, int_x, Two);
        
        signgam = -1;
        for(i = int_x->exp + 1; i < PRECISION; i++)
        {
            if(getMant_UF(int_x, i) != 0)
            {
                signgam = 1;
                break;
            }
        }
        delete_UFs(Two, int_x, new_x, NULL);
        return huge_val;
    }
    
    if(isZero_UF(x))
    {
        signgam = 1;
        delete_UFs(Two, new_x, NULL);
        return huge_val;
    }

    delete_UF(huge_val);
    
    if(isNan_UF(x))
    {
        res = createZero_UF();
        res->kind = NaN;
        delete_UFs(Two, new_x, NULL);
        return res;
    }

    One = createOne_UF();
    if(compare_UF(x, One) == 0 || compare_UF(x, Two) == 0)
    {
        res = createZero_UF();
        res->sign = 1;
        signgam = 1;
        delete_UFs(One, Two, new_x, NULL);
        return res;
    }

    if(isInfinity_UF(x))
    {
        res = createZero_UF();
        res->kind = Infinity;
        res->sign = 1;
        if(x->sign == 1)
            signgam = 1;
        delete_UFs(One, Two, new_x, NULL);
        return res;
    }
    
    //change argument
    cstr = create_CString("10");
    Ten = convertString_UF(cstr);
    delete_CString(cstr);
    z = createZero_UF();
    if(compare_UF(x, Ten) == -1)
    {
        delete_UF(z);
        z = sub_UF(Ten, x);
        if(z->exp>0)
        {
            for(i = z->exp + 1; i < PRECISION + 1; i++)
                setMant_UF(z, i, 0);
        }
        else 
        {
            delete_UF(z);
            z = createZero_UF();
        }
        z = call2_arg1(&add_UF, z, One);
    
        delete_UF(new_x);
        new_x = add_UF(x, z);
        num = convertUnifloat_Integer(z, &error);
    
        delete_UF(z);
        z = clone(x);
        for(i = 1; i < num; i++)
        {
            optmp = convertInteger_UF(i);
            optmp = call2_arg2(&add_UF, x,  optmp);
            z = call2_arg1(&mul_UF, z, optmp);
        }
    }
    delete_UF(Ten);

    Half = createOne_UF();
    Half->exp--;
    res = clone(new_x);
    res = call2_arg1(&sub_UF, res, Half);
    
    optmp = log_UF(new_x);
    res = call2_arg1(&mul_UF, res, optmp);
    delete_UF(optmp);
    res = call2_arg1(&sub_UF, res, new_x);
    optmp = mul_UF(Pi, Two);
    optmp = call1_arg1(&log_UF, optmp);
    optmp = call2_arg1(&mul_UF, optmp, Half);
    res = call2_arg1(&add_UF, res, optmp);
    delete_UF(optmp);
    optmp = gammaSeries_UF(new_x);
    res = call2_arg1(&add_UF, res, optmp);
    delete_UF(optmp);
    
    if(!isZero_UF(z))
    {
        optmp = log_UF(z);
        res = call2_arg1(&sub_UF, res, optmp);
        delete_UF(optmp);
    }

    signgam = 1;
    delete_UFs(Half, Two, One, new_x, NULL);
    return res;
}

Unifloat* lgamma_UF(Unifloat* x)
{
    Unifloat* res, * Half, * One, * Two, * int_x, * Ten, * z, * new_x, * optmp;
    Unifloat* huge_val = clone(infinity_UF);
    int i, num, error;
    CString * cstr;

    Two = createOne_UF();
    Two->exp++;
    new_x = clone(x);

    if((x->sign == -1) && isNormal_UF(x))
    {
        int_x = clone(x);
        for(i = int_x->exp + 1; i < PRECISION; i++)
            setMant_UF(int_x, i, 0);

        int_x = call2_arg1(&div_UF, int_x, Two);
        
        signgam = -1;
        for(i = int_x->exp + 1; i < PRECISION; i++)
        {
            if(getMant_UF(int_x, i) != 0)
            {
                signgam = 1;
                break;
            }
        }

        delete_UFs(Two, int_x, NULL);
        return huge_val;
    }

    if(isZero_UF(x))
    {
        signgam = 1;
        delete_UF(Two);
        return huge_val;
    }
    
    if(isNan_UF(x))
    {
        res = createZero_UF();
        res->kind = NaN;
        delete_UF(Two);
        return res;
    }
    
    One = createOne_UF();
    if(compare_UF(x, One) == 0 || compare_UF(x, Two) == 0)
    {
        res = createZero_UF();
        res->sign = 1;
        signgam = 1;
    
        delete_UFs(One, Two, NULL);
        return res;
    }

    if(isInfinity_UF(x))
    {
        res = createZero_UF();
        res->kind = Infinity;
        res->sign = 1;
        if(x->sign == 1)
            signgam = 1;
        
        delete_UFs(One, Two, NULL);
        return res;
    }
    
    //change argument
    cstr = create_CString("10");
    Ten = convertString_UF(cstr);
    delete_CString(cstr);
    z = createZero_UF();
    if(compare_UF(x, Ten) == -1)
    {
        delete_UF(z);
        z = sub_UF(Ten, x);
        if(z->exp>0)
        {
            for(i = z->exp + 1; i < PRECISION + 1; i++)
                setMant_UF(z, i, 0);
        }
        else 
        {
            delete_UF(z);
            z = createZero_UF();
        }
        z = call2_arg1(&add_UF, z, One);
        delete_UF(new_x);
        new_x = add_UF(x, z);
        num = convertUnifloat_Integer(z, &error);
        delete_UF(z);
        z = clone(x);
        for(i = 1; i < num; i++)
        {
            optmp = convertInteger_UF(i);
            optmp = call2_arg2(&add_UF, x,  optmp);
            z = call2_arg1(&mul_UF, z, optmp);
        }
    }
    delete_UF(Ten);

    Half = createOne_UF();
    Half->exp--;

    res = clone(new_x);
    res = call2_arg1(&sub_UF, res, Half);
    
    optmp = log_UF(new_x);
    res = call2_arg1(&mul_UF, res, optmp);
    delete_UF(optmp);
    
    res = call2_arg1(&sub_UF, res, new_x);
    
    optmp = mul_UF(Pi, Two);
    optmp = call1_arg1(&log_UF, optmp);
    optmp = call2_arg1(&mul_UF, optmp, Half);
    res = call2_arg1(&add_UF, res, optmp);
    delete_UF(optmp);
    
    optmp = gammaSeries_UF(new_x);
    res = call2_arg1(&add_UF, res, optmp);
    delete_UF(optmp);
    
    if(!isZero_UF(z))
    {
        optmp = log_UF(z);
        res = call2_arg1(&sub_UF, res, optmp);
        delete_UF(optmp);
    }
    
    signgam = 1;
    
    delete_UFs(Half, One, Two, new_x, NULL);
    return res;
}

Unifloat* tgamma_UF(Unifloat* x)
{
    Unifloat* res, * Half, * One, * Two, * Ten, * z, * new_x, * optmp;
    Unifloat* huge_val = clone(infinity_UF);
    int i, num, error;

    One = createOne_UF();
    Two = createOne_UF();
    Two->exp++;
    Ten = convertString_UF(create_CString("10"));
    new_x = clone(x);
    z = createZero_UF();
    if((x->sign == -1) && isNormal_UF(x))
    {
        res = createZero_UF();
        res->kind = NaN;
        return res;
    }

    if(isNan_UF(x))
    {
        res = createZero_UF();
        res->kind = NaN;
        return res;
    }

    if(isInfinity_UF(x) && (x->sign == 1))
    {
        res = createZero_UF();
        res->kind = Infinity;
        res->sign = 1;
        return res;
    }

    if(isZero_UF(x))
        return huge_val;

    if(isInfinity_UF(x) && (x->sign == -1))
    {
        res = createZero_UF();
        res->kind = NaN;
        return res;
    }

    //change argument
    if(compare_UF(x, Ten) == -1)
    {
        z = sub_UF(Ten, x);
        if(z->exp>0)
        {
            for(i = z->exp + 1; i < PRECISION + 1; i++)
                setMant_UF(z, i, 0);
        }
        else
            z = createZero_UF();

        z = add_UF(z, One);
        new_x = add_UF(x, z);
        num = convertUnifloat_Integer(z, &error);
        z = clone(x);
        for(i = 1; i < num; i++)
            z = mul_UF(z, add_UF(x, convertInteger_UF(i)));
    }

    Half = createOne_UF();
    Half->exp--;

    res = clone(new_x);
    res = call2_arg1(&sub_UF, res, Half);
    optmp = log_UF(new_x);
    res = mul_UF(res, optmp);
    delete_UF(optmp);
    res = call2_arg1(&sub_UF, res, new_x);
    
    optmp = mul_UF(Pi, Two);
    optmp = call1_arg1(&log_UF, optmp);
    optmp = call2_arg1(&mul_UF, optmp, Half);
    res = call2_arg1(&add_UF, res, optmp);
    delete_UF(optmp);
    
    optmp = gammaSeries_UF(new_x);
    res = call2_arg1(&add_UF, res, optmp);
    
    if(!isZero_UF(z))
    {
        optmp = log_UF(z);
        res = call2_arg1(&sub_UF, res, optmp);
        delete_UF(optmp);
    }
    res = call1_arg1(&exp_UF, res);

    return res;
}

/********************************************************************/
/**                     Unifloat gamma functions                   **/
/********************************************************************/

Unifloat* gammaSeries_UF(Unifloat* x)
{
    CString *cstr;
    Unifloat* result = createZero_UF();
    Unifloat* current1 = createZero_UF();
    Unifloat* current2 = createZero_UF();
    Unifloat* A;
    Unifloat* x_sq;
    x_sq = mul_UF(x, x);

    ////////
    current1 = convertString_UF(create_CString("-261082718496449122051"));
    current2 = convertString_UF(create_CString("21106800"));
    A = div_UF(current1, current2);
    result = div_UF(A, x_sq);

    ////////
    current1 = convertString_UF(create_CString("154210205991661"));
    current2 = convertString_UF(create_CString("444"));
    A = div_UF(current1, current2);
    result = add_UF(A, result);
    result = div_UF(result, x_sq);

    ////////
    cstr = create_CString("-26315271553053477373");
    current1 = convertString_UF(cstr);
    current2 = convertString_UF(create_CString("2418179400"));
    A = div_UF(current1, current2);
    result = add_UF(A, result);
    result = div_UF(result, x_sq);

    ////////
    current1 = convertString_UF(create_CString("151628697551"));
    current2 = convertString_UF(create_CString("396"));
    A = div_UF(current1, current2);
    result = add_UF(A, result);
    result = div_UF(result, x_sq);

    ////////
    current1 = convertString_UF(create_CString("-7709321041217"));
    current2 = convertString_UF(create_CString("505920"));
    A = div_UF(current1, current2);
    result = add_UF(A, result);
    result = div_UF(result, x_sq);

    ////////
    current1 = convertString_UF(create_CString("1723168255201"));
    current2 = convertString_UF(create_CString("2492028"));
    A = div_UF(current1, current2);
    result = add_UF(A, result);
    result = div_UF(result, x_sq);

    ////////
    current1 = convertString_UF(create_CString("-3392780147"));
    current2 = convertString_UF(create_CString("93960"));
    A = div_UF(current1, current2);
    result = add_UF(A, result);
    result = div_UF(result, x_sq);

    ////////
    current1 = convertString_UF(create_CString("657931"));
    current2 = convertString_UF(create_CString("300"));
    A = div_UF(current1, current2);
    result = add_UF(A, result);
    result = div_UF(result, x_sq);

    ////////
    current1 = convertString_UF(create_CString("-236364091"));
    current2 = convertString_UF(create_CString("1506960"));
    A = div_UF(current1, current2);
    result = add_UF(A, result);
    result = div_UF(result, x_sq);

    ////////
    current1 = convertString_UF(create_CString("77683"));
    current2 = convertString_UF(create_CString("5796"));
    A = div_UF(current1, current2);
    result = add_UF(A, result);
    result = div_UF(result, x_sq);

    ////////
    current1 = convertString_UF(create_CString("-174611"));
    current2 = convertString_UF(create_CString("125400"));
    A = div_UF(current1, current2);
    result = add_UF(A, result);
    result = div_UF(result, x_sq);

    ////////
    current1 = convertString_UF(create_CString("43867"));
    current2 = convertString_UF(create_CString("244188"));
    A = div_UF(current1, current2);
    result = add_UF(A, result);
    result = div_UF(result, x_sq);

    ////////
    current1 = convertString_UF(create_CString("-3617"));
    current2 = convertString_UF(create_CString("122400"));
    A = div_UF(current1, current2);
    result = add_UF(A, result);
    result = div_UF(result, x_sq);

    ////////
    current1 = convertString_UF(create_CString("1"));
    current2 = convertString_UF(create_CString("156"));
    A = div_UF(current1, current2);
    result = add_UF(A, result);
    result = div_UF(result, x_sq);

    ////////
    current1 = convertString_UF(create_CString("-691"));
    current2 = convertString_UF(create_CString("360360"));
    A = div_UF(current1, current2);
    result = add_UF(A, result);
    result = div_UF(result, x_sq);

    ////////
    current1 = convertString_UF(create_CString("1"));
    current2 = convertString_UF(create_CString("1188"));
    A = div_UF(current1, current2);
    result = add_UF(A, result);
    result = div_UF(result, x_sq);

    ////////
    current1 = convertString_UF(create_CString("-1"));
    current2 = convertString_UF(create_CString("1680"));
    A = div_UF(current1, current2);
    result = add_UF(A, result);
    result = div_UF(result, x_sq);

    ////////
    current1 = convertString_UF(create_CString("1"));
    current2 = convertString_UF(create_CString("1260"));
    A = div_UF(current1, current2);
    result = add_UF(A, result);
    result = div_UF(result, x_sq);

    ////////
    current1 = convertString_UF(create_CString("-1"));
    current2 = convertString_UF(create_CString("360"));
    A = div_UF(current1, current2);
    result = add_UF(A, result);
    result = div_UF(result, x_sq);

    ////////
    current1 = convertString_UF(create_CString("1"));
    current2 = convertString_UF(create_CString("12"));
    A = div_UF(current1, current2);
    result = add_UF(A, result);
    result = div_UF(result, x);
    ////////

    return result;
}
