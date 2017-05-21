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
#include "unifloat/unifloat_complex.h"
#include "unifloat/bessel.h"
#include "unifloat/exp.h"
#include "unifloat/trig.h"
#include "unifloat/debug.h"

//////////////////////////////////////////////////////////////////////
//                        j0_UF                                     //
//////////////////////////////////////////////////////////////////////
Unifloat* j0_UF(Unifloat* x)
{
    Unifloat* One, * BoundFirst, * BoundSecond, * absx;
    Unifloat* res;

    if(isNan_UF(x))
        return clone(nan_UF);

    if(isInfinity_UF(x))
        return createZero_UF();

    if(isZero_UF(x))
        return createOne_UF();

    One = createOne_UF();
    BoundFirst = clone(One);
    BoundFirst->exp += 3;
    BoundSecond = clone(One);
    BoundSecond->exp += 6;

    absx = abs_UF(x);
    
    if(compare_UF(absx, BoundFirst) != 1)
        res = jnPowerSeries_UF(0, x);
    else if(compare_UF(absx, BoundSecond) != 1)
        res = jnSteed_UF(0, x);
    else
        res = jnHankel_UF(0, x);

    delete_UFs(One, BoundFirst, BoundSecond, absx, NULL);
    return res;
}

//////////////////////////////////////////////////////////////////////
//                        j1_UF                                     //
//////////////////////////////////////////////////////////////////////
Unifloat* j1_UF(Unifloat* x)
{
    Unifloat* One, * BoundFirst, * BoundSecond;
    Unifloat* res = createZero_UF();

    if(isNan_UF(x))
    {
        res->kind = NaN;
        return res;
    }
    
    if(isInfinity_UF(x))
        return res;

    if(isZero_UF(x))
    {
        res->sign = x->sign;
        return res;
    }
    
    One = createOne_UF();
    BoundFirst = clone(One);
    BoundFirst->exp += 3;
    BoundSecond = clone(One);
    BoundSecond->exp += 6;

    if(compare_UF(abs_UF(x), BoundFirst) != 1)
        return jnPowerSeries_UF(1, x);
    else if(compare_UF(abs_UF(x), BoundSecond) != 1)
        return jnSteed_UF(1, x);
    else
        return jnHankel_UF(1, x);
}

//////////////////////////////////////////////////////////////////////
//                        jn_UF                                     //
//////////////////////////////////////////////////////////////////////
Unifloat* jn_UF(int n, Unifloat* x)
{
    int sign = 1, n_sign = 1;
    Unifloat* res = createZero_UF();

    if(n == 0)
        return j0_UF(x);
    else if(n == 1)
        return j1_UF(x);
    else if(n == -1)
    {
        res = j1_UF(x);
        res->sign = (-1)*(res->sign);
        return res;
    }

    if(isNan_UF(x))
    {
        res->kind = NaN;
        return res;
    }
    
    if(isInfinity_UF(x))
        return res;

    if(isZero_UF(x))
    {
        if(((x->sign == 1) && (n > 0)) || ((x->sign == -1) && (n < 0)))
        {
            res->sign = 1;
            return res;
        }
        
        if(((x->sign == -1) && (n > 0)) || ((x->sign == 1) && (n < 0)))
        {
            if(n%2 == 0)
                res->sign = 1;
            else
                res->sign = -1;
            return res;
        }
    }
    //////J[-n]=(-1)^n*J[n]//////
    if(n < 0)
    {
        n_sign = -1;
        if(n % 2 != 0)
            sign = -1;
    }

    if(x->sign == -1)
    {
        if(n % 2 != 0)
            sign = (-1) * sign;
    }
    
    //////compute function//////
    switch(besselMethod_UF(n_sign * n, abs_UF(x)))
    {
    case 1: res = jnPowerSeries_UF(n_sign * n, abs_UF(x));
            break;
    case 2: res = jnHankel_UF(n_sign * n, abs_UF(x));
            break;
    case 3: res = jnMeisselFirst_UF(n_sign * n, abs_UF(x));
            break;
    case 4: res = jnMeisselSecond_UF(n_sign * n, abs_UF(x));
            break;
    case 5: res = jnSteed_UF(n_sign * n, abs_UF(x));
            break;
    case 6: res = jnRecurrent_UF(n_sign * n, abs_UF(x));
            break;
    }
    
    res->sign = sign * res->sign;
    return res;
}

//////////////////////////////////////////////////////////////////////
//                        y0_UF                                     //
//////////////////////////////////////////////////////////////////////
Unifloat* y0_UF(Unifloat* x)
{
    Unifloat* One, * BoundFirst, * BoundSecond;
    Unifloat* res = createZero_UF();

    if(isNan_UF(x))
    {
        res->kind = NaN;
        return res;
    }
    
    if(isInfinity_UF(x) && x->sign == 1)
        return res;

    if(isZero_UF(x) || 
       ((isNormal_UF(x) || isInfinity_UF(x)) && x->sign == -1))
    {
       res->kind = Infinity;
       res->sign = -1;
       return res;
    }

    One = createOne_UF();
    BoundFirst = clone(One);
    BoundFirst->exp += 3;
    BoundSecond = clone(One);
    BoundSecond->exp += 6;

    if(compare_UF(abs_UF(x), BoundFirst) != 1)
        return ynPowerSeries_UF(0, x);
    else if(compare_UF(abs_UF(x), BoundSecond) != 1)
        return ynSteed_UF(0, x);
    else
        return ynHankel_UF(0, x);
}


//////////////////////////////////////////////////////////////////////
//                        y1_UF                                     //
//////////////////////////////////////////////////////////////////////
Unifloat* y1_UF(Unifloat* x)
{
    Unifloat* One, * BoundFirst, * BoundSecond;
    Unifloat* res = createZero_UF();

    if(isNan_UF(x))
    {
        res->kind = NaN;
        return res;
    }
    
    if(isInfinity_UF(x) && x->sign == 1)
        return res;

    if(isZero_UF(x) || 
       ((isNormal_UF(x) || isInfinity_UF(x)) && x->sign == -1))
    {
       res->kind = Infinity;
       res->sign = -1;
       return res;
    }
    
    One = createOne_UF();
    BoundFirst = clone(One);
    BoundFirst->exp += 3;
    BoundSecond = clone(One);
    BoundSecond->exp += 6;

    if(compare_UF(abs_UF(x), BoundFirst) != 1)
        return ynPowerSeries_UF(1, x);
    else if(compare_UF(abs_UF(x), BoundSecond) != 1)
        return ynSteed_UF(1, x);
    else
        return ynHankel_UF(1, x);
}

//////////////////////////////////////////////////////////////////////
//                        yn_UF                                     //
//////////////////////////////////////////////////////////////////////
Unifloat* yn_UF (int n, Unifloat* x)
{
    int sign = 1, n_sign = 1;
    Unifloat* res = createZero_UF();

    if(n == 0)
        return y0_UF(x);
    if((n == 1) || (n == -1))
        return y1_UF(x);

    if(isNan_UF(x))
    {
        res->kind = NaN;
        return res;
    }
    
    if(isInfinity_UF(x) && x->sign == 1)
        return res;

    if(isZero_UF(x) || 
       ((isNormal_UF(x) || isInfinity_UF(x)) && x->sign == -1))
    {
       res->kind = Infinity;
       res->sign = -1;
       return res;
    }
    
    //////Y[-n]=(-1)^n*Y[n]//////
    if(n < 0)
    {
        n_sign = -1;
        if(n % 2 != 0)
            sign = -1;
    }

    //////compute function//////
    switch(besselMethod_UF(n_sign * n, abs_UF(x)))
    {
    case 1: res = ynPowerSeries_UF(n_sign * n, x);
            break;
    case 2: res = ynHankel_UF(n_sign * n, x);
            break;
    case 3: res = ynMeisselFirst_UF(n_sign * n, x);
            break;
    case 4: res = ynMeisselSecond_UF(n_sign * n, x);
            break;
    case 5: res = ynSteed_UF(n_sign * n, x);
            break;
    case 6: res = ynRecurrent_UF(n_sign * n, abs_UF(x));
            break;
    }

    res->sign = sign * res->sign;
    return res;
}

//////////////////////////////////////////////////////////////////////
//                    jnPowerSeries_UF                              //
//////////////////////////////////////////////////////////////////////
Unifloat* jnPowerSeries_UF(int n, Unifloat* x)
{
    ////// calculate through representation function in power sum //////
    int num;
    Unifloat* res, * z, * One;
    Unifloat* x_pow_n;
    Unifloat* n_fact;
    Unifloat* tmp1, * tmp2, * tmp3, * tmp;

    One = createOne_UF();
    z = mul_UF(x, x);
    res = createOne_UF();

    for(num = 45; num > 0; num--)
    {
        tmp1 = mul_UF(res, z);
        delete_UF(res);

        tmp2 = convertInteger_UF(num * (n + num));
        tmp3 = div_UF(tmp1, tmp2);
        tmp3->exp -= 2;

        res = sub_UF(One, tmp3);

        delete_UFs(tmp1, tmp2, tmp3, NULL);
    }

    if((n != 0) && (n <= 100))
    {
        x_pow_n = power_UF(x, n);
        n_fact = factorial_UF(n);

        tmp = mul_UF(res, x_pow_n);
        delete_UF(res);

        res = div_UF(tmp, n_fact);
        res->exp -= n;

        delete_UFs(x_pow_n, n_fact, tmp, NULL);
    }

    delete_UFs(z, One, NULL);

    return res;
}

//////////////////////////////////////////////////////////////////////
//                    ynPowerSeries_UF                              //
//////////////////////////////////////////////////////////////////////
Unifloat* ynPowerSeries_UF(int n, Unifloat* x)
{
    int num;
    Unifloat* currentsumA, * aux, * z;
    Unifloat* currentsumB = createZero_UF();
    Unifloat* work = createZero_UF();
    Unifloat* currentA, * currentB;
    Unifloat* M = createOne_UF();

    z = mul_UF(x, x);

    ///////////////////////////////////////////////////////////////////////////
    //first part of representation
    ///////////////////////////////////////////////////////////////////////////
    aux = clone(x);
    aux->exp--;
    currentsumA = add_UF(log_UF(aux), Gamma);
    currentsumA = mul_UF(jnPowerSeries_UF(n, x), currentsumA);
    currentsumA->exp++;
        
    ///////////////////////////////////////////////////////////////////////////
    //second part of representation
    ///////////////////////////////////////////////////////////////////////////
    if(n != 0)
    {
        currentB = createOne_UF();
        currentsumB = createOne_UF();
        for(num = 1; num < n; num++)
        {
            currentB = mul_UF(currentB, z);
            currentB = div_UF(currentB, convertInteger_UF(num*(n - num)));
            currentB->exp -= 2;
            currentsumB = add_UF(currentsumB, currentB);
        }
        
        aux = createOne_UF();
        for(num = 1; num < n; num++)
        {
            M = mul_UF(M, convertInteger_UF(num));
            aux = mul_UF(aux, x);
        }
        aux = mul_UF(aux, x);
        M = div_UF(M, aux);
        currentsumB = mul_UF(currentsumB, M);
        currentsumB->exp += n;

        currentsumA = sub_UF(currentsumA, currentsumB);
    }

    ///////////////////////////////////////////////////////////////////////////
    //third part of representation
    ///////////////////////////////////////////////////////////////////////////
    currentsumB = createZero_UF();
    currentB = createOne_UF();
    currentA = createOne_UF();
    for(num = 1; num < n + 1; num++)
        work = add_UF(work, div_UF(createOne_UF(), convertInteger_UF(num)));

    currentsumB = add_UF(currentsumB, work);
    num = 1;
    while(currentB->exp > -PRECISION - 1)
    {
        work = add_UF(work, div_UF(createOne_UF(), convertInteger_UF(num)));
        work = add_UF(work, div_UF(createOne_UF(), convertInteger_UF(n + num)));
        currentA = mul_UF(currentA, z);
        currentA = div_UF(currentA, convertInteger_UF(num * (n + num)));
        currentA->exp -= 2;
            
        currentA->sign = (-1)*currentA->sign;
        currentB = mul_UF(currentA, work);
        currentsumB = add_UF(currentsumB, currentB);
        num++;
    }

    if(n != 0)
    {
        M = mul_UF(M, convertInteger_UF(n));
        currentsumB->exp -= n;
        currentsumB = div_UF(currentsumB, M);
    }

    currentsumA = sub_UF(currentsumA, currentsumB);
    currentsumA = div_UF(currentsumA, Pi);
    
    return currentsumA;
}

//////////////////////////////////////////////////////////////////////
//                        jnHankel_UF                               //
//////////////////////////////////////////////////////////////////////
Unifloat* jnHankel_UF(int n, Unifloat* x)
{
    ////// asimptotic representation //////
    int num = 1, exp1 = 1, exp2 = 0;   

    Unifloat* Four, * INDEX, * z, * aux2;
    Unifloat* aux1 = createZero_UF();
    Unifloat* currentsumA = createZero_UF();
    Unifloat* currentsumB = createZero_UF();
    Unifloat* currentA, * currentB, * arg, * pi_4;
    Unifloat* One, * Two;

    One = createOne_UF();
    Two = createOne_UF();
    Two->exp++;
    Four = createOne_UF();
    Four->exp += 2;
    INDEX = convertInteger_UF(n);
    pi_4 = clone(Pi);
    pi_4->exp -= 2;
    z = mul_UF(x, x);
    z->exp += 6;
    
    //////calculating argument of sin() and cos()//////
    arg = mul_UF(Pi, INDEX);
    arg->exp--;
    arg = sub_UF(x, arg);
    arg = sub_UF(arg, pi_4);

    //////calculating all//////
    currentA = createOne_UF();
    aux2 = mul_UF(INDEX,INDEX);
    currentB = clone(aux2);
    aux2->exp += 2;
    currentB->exp += 2;
    currentB = sub_UF(currentB, One);
    currentB = div_UF(currentB, x);
    currentB->exp -= 3;
    
    while((exp2 <= exp1) && (exp2 > -PRECISION - 1))
    {
        if(currentB->exp < 0)
            exp1 = currentB->exp;
        
        currentsumA = add_UF(currentsumA, currentA);
        currentsumB = add_UF(currentsumB, currentB);

        aux1 = clone(aux2);
        aux1 = sub_UF(aux1, convertInteger_UF((4*num-3)*(4*num-3)));
        currentA = mul_UF(currentA, aux1);
        aux1 = clone(aux2);
        aux1 = sub_UF(aux1, convertInteger_UF((4*num-1)*(4*num-1)));
        currentA = mul_UF(currentA, aux1);
        aux1 = convertInteger_UF(2*num*(2*num-1));
        aux1 = mul_UF(aux1, z);
        currentA = div_UF(currentA, aux1);

        aux1 = clone(aux2);
        aux1 = sub_UF(aux1, convertInteger_UF((4*num-1)*(4*num-1)));
        currentB = mul_UF(currentB, aux1);
        
        aux1 = clone(aux2);
        aux1 = sub_UF(aux1, convertInteger_UF((4*num+1)*(4*num+1)));
        currentB = mul_UF(currentB, aux1);
        
        aux1 = convertInteger_UF(2*num*(2*num+1));
        aux1 = mul_UF(aux1, z);
        currentB = div_UF(currentB, aux1);

        currentA->sign=(-1)*currentA->sign;
        currentB->sign=(-1)*currentB->sign;
        
        if(currentB->exp < 0)
            exp2 = currentB->exp;
        
        num++;
    }
    
    currentsumA = mul_UF(currentsumA, cos_UF(arg));
    currentsumB = mul_UF(currentsumB, sin_UF(arg));
    currentsumA = sub_UF(currentsumA, currentsumB);
    currentsumA = mul_UF(currentsumA, sqrt_UF(div_UF(Two, mul_UF(Pi, x))));
    
    return currentsumA;
}

//////////////////////////////////////////////////////////////////////
//                        ynHankel_UF                               //
//////////////////////////////////////////////////////////////////////
Unifloat* ynHankel_UF(int n, Unifloat* x)
{
    ////// asimptotic representation //////
    int num = 1, exp1 = 1, exp2 = 0;

    Unifloat* Four, * INDEX, * z, * aux2;
    Unifloat* aux1 = createZero_UF();
    Unifloat* currentsumA = createZero_UF();
    Unifloat* currentsumB = createZero_UF();
    Unifloat* currentA, * currentB, * arg;

    Four = convertInteger_UF(4);
    INDEX = convertInteger_UF(n);
    z = mul_UF(x, mul_UF(x, convertInteger_UF(64)));
    
    //////calculating argument of sin() and cos()//////
    arg = sub_UF(x, mul_UF(div_UF(Pi, convertInteger_UF(2)), INDEX));
    arg = sub_UF(arg, div_UF(Pi, Four));
    aux2 = mul_UF(INDEX,INDEX);
    aux2 = mul_UF(aux2,Four);

    //////calculating all//////
    currentA = convertInteger_UF(1);
    currentB = mul_UF(INDEX, INDEX);
    currentB = mul_UF(currentB, Four);
    currentB = sub_UF(currentB, convertInteger_UF(1));
    currentB = div_UF(currentB, x);
    currentB = div_UF(currentB, convertInteger_UF(8));
    while((exp2 <= exp1)&&(exp2>-PRECISION-1))
    {
        if(currentB->exp<0)
            exp1 = currentB->exp;
        currentsumA = add_UF(currentsumA, currentA);
        currentsumB = add_UF(currentsumB, currentB);

        aux1 = clone(aux2);
        aux1 = sub_UF(aux1, mul_UF(convertInteger_UF(4*num-3),convertInteger_UF(4*num-3)));
        currentA = mul_UF(currentA, aux1);
        aux1 = clone(aux2);
        aux1 = sub_UF(aux1, mul_UF(convertInteger_UF(4*num-1),convertInteger_UF(4*num-1)));
        currentA = mul_UF(currentA, aux1);
        aux1 = mul_UF(convertInteger_UF(2*num), convertInteger_UF(2*num-1));
        aux1 = mul_UF(aux1, z);
        currentA = div_UF(currentA, aux1);

        aux1 = clone(aux2);
        aux1 = sub_UF(aux1, mul_UF(convertInteger_UF(4*num-1),convertInteger_UF(4*num-1)));
        currentB = mul_UF(currentB, aux1);
        aux1 = clone(aux2);
        aux1 = sub_UF(aux1, mul_UF(convertInteger_UF(4*num+1),convertInteger_UF(4*num+1)));
        currentB = mul_UF(currentB, aux1);
        aux1 = mul_UF(convertInteger_UF(2*num), convertInteger_UF(2*num+1));
        aux1 = mul_UF(aux1, z);
        currentB = div_UF(currentB, aux1);

        currentA->sign=(-1)*currentA->sign;
        currentB->sign=(-1)*currentB->sign;
        if(currentB->exp<0)
            exp2 = currentB->exp;
        
        num++;
    }
    
    currentsumA = mul_UF(currentsumA, sin_UF(arg));
    currentsumB = mul_UF(currentsumB, cos_UF(arg));
    currentsumA = add_UF(currentsumA, currentsumB);
    currentsumA = mul_UF(currentsumA, sqrt_UF(div_UF(convertInteger_UF(2), mul_UF(Pi, x))));
    
    return currentsumA;
}

//////////////////////////////////////////////////////////////////////
//                          jnSteed_UF                              //
//////////////////////////////////////////////////////////////////////
Unifloat* jnSteed_UF(int n, Unifloat* x)
{
    int num, N1 = 1, N2 = 300, sign;
    Unifloat* Wronskian, * J, * gamma, * CF1, * A1, * A2, * B1, * B2, * AA;
    Unifloat* One, * Two, * Four, * INDEX, * Zero;
    Unifloat* temp, * tmp, * Q0, * Q1;
    UnifloatComplex* CF2, * aux;

    Zero = createZero_UF();
    One = createOne_UF();
    Two = createOne_UF();
    Two->exp++;
    Four = createOne_UF();
    Four->exp += 2;
    INDEX = convertInteger_UF(n);

    CF2 = create_UFComplex(clone(Zero), clone(Zero));
    aux = create_UFComplex(clone(Zero), clone(Zero));

    //////determine N1///////
    B1 = div_UF(createOne_UF(), x);
    B1->exp++;
    A1 = mul_UF(INDEX, B1);
    Q0 = clone(A1);
    Q1 = mul_UF(Q0, add_UF(A1, B1));
    Q1 = sub_UF(Q1, One);
    temp = add_UF(A1, B1);
    while(Q1->exp < 50)
    {
        N1 += 1;
        temp = add_UF(temp, B1);
        tmp = mul_UF(temp, Q1);
        tmp = sub_UF(tmp, Q0);
        Q0 = clone(Q1);
        Q1 = clone(tmp);
    }
    
    //////compute CF1//////
    A1 = clone(One);
    A2 = clone(Zero);
    
    B1 = clone(Zero);
    B2 = clone(One);

    for(num = 1; num < N1; num++)
    {
        AA = clone(A2);
        A2 = sub_UF(div_UF(mul_UF(mul_UF(Two, A2), 
           add_UF(INDEX, convertInteger_UF(num))), x), A1);
        A1 = clone(AA);
        AA = clone(B2);
        B2 = sub_UF(div_UF(mul_UF(mul_UF(Two, B2), 
           add_UF(INDEX, convertInteger_UF(num))), x), B1);
        B1 = clone(AA);
    }
    
    CF1 = div_UF(A2, B2);
    CF1 = add_UF(div_UF(INDEX, x), CF1);
    sign = B2->sign;

    CF2->Re = clone(x);
    CF2->Re->exp++;
    CF2->Im = convertInteger_UF(N2 + 1);
    CF2->Im->exp++;
    for(num = N2; num > 0; num--)
    {
        aux->Re = mul_UF(convertInteger_UF(2*num+1), convertInteger_UF(2*num+1));
        aux->Re->exp -= 2;
        aux->Re = sub_UF(aux->Re, mul_UF(INDEX, INDEX));
        aux->Im = createZero_UF();
        CF2 = div_UFComplex(aux, CF2);
        aux->Re = clone(x);
        aux->Re->exp++;
        aux->Im = convertInteger_UF(num);
        aux->Im->exp++;
        CF2 = add_UFComplex(aux, CF2);
    }
    aux->Im = clone(One);
    aux->Im->exp -= 2;
    aux->Im = sub_UF(aux->Im, mul_UF(INDEX, INDEX));
    aux->Im = div_UF(aux->Im, x);
    aux->Re = createZero_UF();

    CF2 = div_UFComplex(aux, CF2);
    CF2->Im = add_UF(CF2->Im, One);
    CF2->Re = sub_UF(CF2->Re, div_UF(One, mul_UF(Two, x)));

    //////compute J//////
    Wronskian = div_UF(Two, mul_UF(Pi, x));

    gamma = div_UF(sub_UF(CF2->Re, CF1), CF2->Im);
    J = div_UF(Wronskian, add_UF(CF2->Im, mul_UF(gamma, sub_UF(CF2->Re, CF1))));
    J = sqrt_UF(J);
    J->sign = sign;
    return J;
}

//////////////////////////////////////////////////////////////////////
//                          ynSteed_UF                              //
//////////////////////////////////////////////////////////////////////
Unifloat* ynSteed_UF(int n, Unifloat* x)
{
    int num, N1 = 1, N2 = 300, sign;
    Unifloat* Wronskian, * Y, * gamma, * CF1, * A1, * A2, * B1, * B2, * AA;
    Unifloat* One, * Two, * Four, * INDEX, * Zero;
    Unifloat* temp, * tmp, * Q0, * Q1;
    UnifloatComplex* CF2, * aux;

    Zero = createZero_UF();
    One = convertInteger_UF(1);
    Two = convertInteger_UF(2);
    Four = convertInteger_UF(4);
    INDEX = convertInteger_UF(n);

    CF2 = create_UFComplex(clone(Zero), clone(Zero));
    aux = create_UFComplex(clone(Zero), clone(Zero));

    //////determine N1///////
    B1 = div_UF(createOne_UF(), x);
    B1->exp++;
    A1 = mul_UF(INDEX, B1);
    Q0 = clone(A1);
    Q1 = mul_UF(Q0, add_UF(A1, B1));
    Q1 = sub_UF(Q1, One);
    temp = add_UF(A1, B1);
    while(Q1->exp < 50)
    {
        N1 += 1;
        temp = add_UF(temp, B1);
        tmp = mul_UF(temp, Q1);
        tmp = sub_UF(tmp, Q0);
        Q0 = clone(Q1);
        Q1 = clone(tmp);
    }

    //////compute CF1//////
    A1 = clone(One);
    A2 = clone(Zero);
    
    B1 = clone(Zero);
    B2 = clone(One);

    for(num = 1; num < N1; num++)
    {
        AA = clone(A2);
        A2 = sub_UF(div_UF(mul_UF(mul_UF(Two, A2), 
           add_UF(INDEX, convertInteger_UF(num))), x), A1);
        A1 = clone(AA);
        AA = clone(B2);
        B2 = sub_UF(div_UF(mul_UF(mul_UF(Two, B2), 
           add_UF(INDEX, convertInteger_UF(num))), x), B1);
        B1 = clone(AA);
    }
    
    CF1=div_UF(A2, B2);
    CF1 = add_UF(div_UF(INDEX, x), CF1);
    sign = B2->sign;

    CF2->Re = mul_UF(Two, x);
    CF2->Im = mul_UF(Two, convertInteger_UF(N2 + 1));

    for(num = N2; num > 0; num--)
    {
        aux->Re = div_UF(mul_UF(convertInteger_UF(2*num+1), 
                               convertInteger_UF(2*num+1)), Four);
        aux->Re = sub_UF(aux->Re, mul_UF(INDEX, INDEX));
        aux->Im = convertInteger_UF(0);
        CF2 = div_UFComplex(aux, CF2);
        aux->Re = mul_UF(Two, x);
        aux->Im = mul_UF(Two, convertInteger_UF(num));
        CF2 = add_UFComplex(aux, CF2);
    }
    aux->Im = div_UF(One, Four);
    aux->Im = sub_UF(aux->Im, mul_UF(INDEX, INDEX));
    aux->Im = div_UF(aux->Im, x);
    aux->Re = convertInteger_UF(0);

    CF2 = div_UFComplex(aux, CF2);
    CF2->Im = add_UF(CF2->Im, One);
    CF2->Re = sub_UF(CF2->Re, div_UF(One, mul_UF(Two, x)));

    //////compute Y//////
    Wronskian = div_UF(convertInteger_UF(2), mul_UF(Pi, x));

    gamma = div_UF(sub_UF(CF2->Re, CF1), CF2->Im);
    Y = div_UF(Wronskian, add_UF(CF2->Im, mul_UF(gamma, sub_UF(CF2->Re, CF1))));
    Y = sqrt_UF(Y);
    Y->sign=sign;
    Y = mul_UF(Y, gamma);
    return Y;
}

//////////////////////////////////////////////////////////////////////
//                          besselMethod_UF                         //
//////////////////////////////////////////////////////////////////////
int besselMethod_UF(int n, Unifloat* x)
{
    Unifloat* tmp, * INDEX = convertInteger_UF(n);
    Unifloat* Ten = convertInteger_UF(10);
    Unifloat* One, * BoundFirst, * BoundSecond;
    
    One = createOne_UF();
    BoundFirst = clone(One);
    BoundFirst->exp += 3;
    BoundSecond = clone(One);
    BoundSecond->exp += 6;

    tmp = createOne_UF();
    tmp->exp += 4;
    BoundFirst = add_UF(BoundFirst, tmp);

    if(n <= 10)
    {
        if(compare_UF(abs_UF(x), BoundFirst) != 1)
            return 1;
        else if(compare_UF(abs_UF(x), BoundSecond) != 1)
            return 5;
        else
            return 2;
    }
    else if(n <= 100000)
        return 6; //Recurrence, forward and backward
    else
    {
        tmp = mul_UF(INDEX, INDEX);
        if(compare_UF(x, tmp) != -1)
            return 2; //Hankel's

        tmp = mul_UF(convertString_UF(create_CString("0.98")), INDEX);
        if(compare_UF(x, tmp) != 1)
            return 3; //Meissel's First
        
        tmp = mul_UF(convertString_UF(create_CString("1.03")), INDEX);
        if(compare_UF(x, tmp) != -1)
            return 4; //Meissel's Second
        
        if((compare_UF(x, add_UF(INDEX, Ten)) != 1)  
            && (compare_UF(x, sub_UF(INDEX, Ten)) != -1))
            return 5; //Steed's
    }
    return 0;
}

//////////////////////////////////////////////////////////////////////
//                        factorialSeries_UF                        //
//////////////////////////////////////////////////////////////////////
static Unifloat* factorialSeries_UF(int n)
{
    Unifloat* result = createZero_UF();
    Unifloat* current = createZero_UF();
    Unifloat* index = convertInteger_UF(n);
    Unifloat* square_index = mul_UF(index, index);

    current = div_UF(convertInteger_UF(1), mul_UF(convertInteger_UF(12), index));
    result = add_UF(result, current);

    index = mul_UF(index, square_index);

    current = div_UF(convertInteger_UF(1), mul_UF(convertInteger_UF(360), index));
    result = sub_UF(result, current);

    index = mul_UF(index, square_index);

    current = div_UF(convertInteger_UF(1), mul_UF(convertInteger_UF(1260), index));
    result = add_UF(result, current);

    index = mul_UF(index, square_index);

    current = div_UF(convertInteger_UF(1), mul_UF(convertInteger_UF(1680), index));
    result = sub_UF(result, current);

    index = mul_UF(index, square_index);

    current = div_UF(convertInteger_UF(5), mul_UF(convertInteger_UF(5940), index));
    result = add_UF(result, current);

    index = mul_UF(index, square_index);

    current = div_UF(convertInteger_UF(691), mul_UF(convertInteger_UF(360360), index));
    result = sub_UF(result, current);

    index = mul_UF(index, square_index);

    current = div_UF(convertInteger_UF(7), mul_UF(convertInteger_UF(1092), index));
    result = add_UF(result, current);

    index = mul_UF(index, square_index);

    current = div_UF(convertInteger_UF(3617), mul_UF(convertInteger_UF(122400), index));
    result = sub_UF(result, current);

    index = mul_UF(index, square_index);

    current = div_UF(convertInteger_UF(43867), mul_UF(convertInteger_UF(244188), index));
    result = add_UF(result, current);

    index = mul_UF(index, square_index);

    current = div_UF(convertInteger_UF(174611), mul_UF(convertInteger_UF(125400), index));
    result = sub_UF(result, current);

    return result;
}

//////////////////////////////////////////////////////////////////////
//                        jnMeisselFirst_UF                         //
//////////////////////////////////////////////////////////////////////
Unifloat* jnMeisselFirst_UF(int n, Unifloat* x)
{
    Unifloat* result, * currentV, * current1;
    Unifloat* current2 = createZero_UF();
    Unifloat* sumV = createZero_UF();
    Unifloat* arg, * aux1, * aux2;

    arg = div_UF(x, convertInteger_UF(n));
    aux1 = sqrt_UF(sub_UF(convertInteger_UF(1), mul_UF(arg, arg)));
    aux2 = sub_UF(convertInteger_UF(1), mul_UF(arg, arg));

    //////compute Vsum
    //V1
    currentV = convertInteger_UF(2);
    currentV = add_UF(currentV, mul_UF(convertInteger_UF(3), power_UF(arg, 2)));
    currentV = div_UF(currentV, aux1);
    currentV = div_UF(currentV, aux2);
    currentV = sub_UF(currentV, convertInteger_UF(2));
    currentV = div_UF(currentV, mul_UF(convertInteger_UF(24), convertInteger_UF(n)));
    sumV = add_UF(sumV, currentV);

    //V2
    currentV = mul_UF(convertInteger_UF(4), power_UF(arg, 2));
    currentV = add_UF(currentV, power_UF(arg, 4));
    currentV = div_UF(currentV, power_UF(convertInteger_UF(n), 2));
    currentV = div_UF(currentV, convertInteger_UF(16));
    currentV = div_UF(currentV, power_UF(aux2, 3));
    sumV = sub_UF(sumV, currentV);

    //V3
    currentV = convertInteger_UF(16);
    currentV = sub_UF(currentV, mul_UF(convertInteger_UF(1512), power_UF(arg, 2)));
    currentV = sub_UF(currentV, mul_UF(convertInteger_UF(3654), power_UF(arg, 4)));
    currentV = sub_UF(currentV, mul_UF(convertInteger_UF(375), power_UF(arg, 6)));
    currentV = div_UF(currentV, power_UF(aux2, 4));
    currentV = div_UF(currentV, aux1);
    currentV = sub_UF(currentV, convertInteger_UF(16));
    currentV = div_UF(currentV, mul_UF(convertInteger_UF(5760), power_UF(convertInteger_UF(n), 3)));
    sumV = sub_UF(sumV, currentV);

    //V4
    currentV = mul_UF(convertInteger_UF(32), power_UF(arg, 2));
    currentV = add_UF(currentV, mul_UF(convertInteger_UF(288), power_UF(arg, 4)));
    currentV = add_UF(currentV, mul_UF(convertInteger_UF(232), power_UF(arg, 6)));
    currentV = add_UF(currentV, mul_UF(convertInteger_UF(13), power_UF(arg, 8)));
    currentV = div_UF(currentV, mul_UF(convertInteger_UF(128), power_UF(convertInteger_UF(n), 4)));
    currentV = div_UF(currentV, power_UF(aux2, 6));
    sumV = sub_UF(sumV, currentV);

    //V5
    currentV = mul_UF(convertInteger_UF(67599), power_UF(arg, 10));
    currentV = add_UF(currentV, mul_UF(convertInteger_UF(1914210), power_UF(arg, 8)));
    currentV = add_UF(currentV, mul_UF(convertInteger_UF(4744640), power_UF(arg, 6)));
    currentV = add_UF(currentV, mul_UF(convertInteger_UF(1891200), power_UF(arg, 4)));
    currentV = add_UF(currentV, mul_UF(convertInteger_UF(78720), power_UF(arg, 2)));
    currentV = add_UF(currentV, convertInteger_UF(256));
    currentV = div_UF(currentV, mul_UF(convertInteger_UF(322560), power_UF(convertInteger_UF(n), 5)));
    currentV = div_UF(currentV, power_UF(aux2, 7));
    currentV = div_UF(currentV, aux1);
    currentV = sub_UF(currentV, div_UF(convertInteger_UF(1), mul_UF(convertInteger_UF(1260), power_UF(convertInteger_UF(n), 5))));
    sumV = add_UF(sumV, currentV);

    //V6
    currentV = convertInteger_UF(48);
    currentV = add_UF(currentV, mul_UF(convertInteger_UF(2580), power_UF(arg, 2)));
    currentV = add_UF(currentV, mul_UF(convertInteger_UF(14884), power_UF(arg, 4)));
    currentV = add_UF(currentV, mul_UF(convertInteger_UF(17493), power_UF(arg, 6)));
    currentV = add_UF(currentV, mul_UF(convertInteger_UF(4242), power_UF(arg, 8)));
    currentV = add_UF(currentV, mul_UF(convertInteger_UF(103), power_UF(arg, 10)));
    currentV = mul_UF(currentV, power_UF(arg, 2));
    currentV = div_UF(currentV, mul_UF(convertInteger_UF(192), power_UF(convertInteger_UF(n), 6)));
    currentV = div_UF(currentV, power_UF(aux2, 9));
    sumV = sub_UF(sumV, currentV);

    //V7
    currentV = mul_UF(convertInteger_UF(881664), power_UF(arg, 2));
    currentV = add_UF(currentV, mul_UF(convertInteger_UF(99783936), power_UF(arg, 4)));
    currentV = add_UF(currentV, mul_UF(convertInteger_UF(1135145088), power_UF(arg, 6)));
    currentV = add_UF(currentV, mul_UF(convertInteger_UF(2000000000), power_UF(arg, 8)));
    currentV = add_UF(currentV, mul_UF(convertInteger_UF(884531440), power_UF(arg, 8)));
    currentV = add_UF(currentV, mul_UF(convertInteger_UF(1965889800), power_UF(arg, 10)));
    currentV = add_UF(currentV, mul_UF(convertInteger_UF(318291750), power_UF(arg, 12)));
    currentV = add_UF(currentV, mul_UF(convertInteger_UF(5635995), power_UF(arg, 14)));
    currentV = sub_UF(currentV, convertInteger_UF(2048));
    currentV = div_UF(currentV, mul_UF(convertInteger_UF(3440640), power_UF(convertInteger_UF(n), 7)));
    currentV = div_UF(currentV, power_UF(aux2, 10));
    currentV = div_UF(currentV, aux1);
    currentV = add_UF(currentV, div_UF(convertInteger_UF(1), mul_UF(convertInteger_UF(1680), power_UF(convertInteger_UF(n), 7))));
    sumV = add_UF(sumV, currentV);

    //V8
    currentV = convertInteger_UF(1024);
    currentV = add_UF(currentV, mul_UF(convertInteger_UF(248320), power_UF(arg, 2)));
    currentV = add_UF(currentV, mul_UF(convertInteger_UF(5095936), power_UF(arg, 4)));
    currentV = add_UF(currentV, mul_UF(convertInteger_UF(24059968), power_UF(arg, 6)));
    currentV = add_UF(currentV, mul_UF(convertInteger_UF(34280896), power_UF(arg, 8)));
    currentV = add_UF(currentV, mul_UF(convertInteger_UF(15252048), power_UF(arg, 10)));
    currentV = add_UF(currentV, mul_UF(convertInteger_UF(1765936), power_UF(arg, 12)));
    currentV = add_UF(currentV, mul_UF(convertInteger_UF(23797), power_UF(arg, 14)));
    currentV = mul_UF(currentV, power_UF(arg, 2));
    currentV = div_UF(currentV, mul_UF(convertInteger_UF(4096), power_UF(convertInteger_UF(n), 8)));
    currentV = div_UF(currentV, power_UF(aux2, 12));
    sumV = sub_UF(sumV, currentV);

    //////compute result
    result = div_UF(convertInteger_UF(1), sqrt_UF(aux1));

    current1 = log_UF(div_UF(arg, add_UF(convertInteger_UF(1), aux1)));
    current1 = add_UF(current1, aux1);
    current1 = mul_UF(current1, convertInteger_UF(n));

    current2 = add_UF(current2, sumV);
    current2 = add_UF(current2, factorialSeries_UF(n));

    result = mul_UF(result, exp_UF(sub_UF(current1, current2)));
    result = div_UF(result, sqrt_UF(mul_UF(mul_UF(convertInteger_UF(2), convertInteger_UF(n)), Pi)));
    return result;
}

//////////////////////////////////////////////////////////////////////
//                        ynMeisselFirst_UF                         //
//////////////////////////////////////////////////////////////////////
Unifloat* ynMeisselFirst_UF(int n, Unifloat* x)
{
    Unifloat* result, * currentV, * current1;
    Unifloat* current2 = createZero_UF();
    Unifloat* sumV = createZero_UF();
    Unifloat* arg, * aux1, * aux2;

    arg = div_UF(x, convertInteger_UF(n));
    aux1 = sqrt_UF(sub_UF(convertInteger_UF(1), mul_UF(arg, arg)));
    aux2 = sub_UF(convertInteger_UF(1), mul_UF(arg, arg));
    
    //////compute Vsum
    //V1
    currentV = convertInteger_UF(2);
    currentV = add_UF(currentV, mul_UF(convertInteger_UF(3), power_UF(arg, 2)));
    currentV = div_UF(currentV, aux1);
    currentV = div_UF(currentV, aux2);
    currentV = sub_UF(currentV, convertInteger_UF(2));
    currentV = div_UF(currentV, mul_UF(convertInteger_UF(24), convertInteger_UF(n)));
    sumV = add_UF(sumV, currentV);

    //V2
    currentV = mul_UF(convertInteger_UF(4), power_UF(arg, 2));
    currentV = add_UF(currentV, power_UF(arg, 4));
    currentV = div_UF(currentV, power_UF(convertInteger_UF(n), 2));
    currentV = div_UF(currentV, convertInteger_UF(16));
    currentV = div_UF(currentV, power_UF(aux2, 3));
    sumV = add_UF(sumV, currentV);

    //V3
    currentV = convertInteger_UF(16);
    currentV = sub_UF(currentV, mul_UF(convertInteger_UF(1512), power_UF(arg, 2)));
    currentV = sub_UF(currentV, mul_UF(convertInteger_UF(3654), power_UF(arg, 4)));
    currentV = sub_UF(currentV, mul_UF(convertInteger_UF(375), power_UF(arg, 6)));
    currentV = div_UF(currentV, power_UF(aux2, 4));
    currentV = div_UF(currentV, aux1);
    currentV = sub_UF(currentV, convertInteger_UF(16));
    currentV = div_UF(currentV, mul_UF(convertInteger_UF(5760), power_UF(convertInteger_UF(n), 3)));
    sumV = sub_UF(sumV, currentV);

    //V4
    currentV = mul_UF(convertInteger_UF(32), power_UF(arg, 2));
    currentV = add_UF(currentV, mul_UF(convertInteger_UF(288), power_UF(arg, 4)));
    currentV = add_UF(currentV, mul_UF(convertInteger_UF(232), power_UF(arg, 6)));
    currentV = add_UF(currentV, mul_UF(convertInteger_UF(13), power_UF(arg, 8)));
    currentV = div_UF(currentV, mul_UF(convertInteger_UF(128), power_UF(convertInteger_UF(n), 4)));
    currentV = div_UF(currentV, power_UF(aux2, 6));
    sumV = add_UF(sumV, currentV);

    //V5
    currentV = mul_UF(convertInteger_UF(67599), power_UF(arg, 10));
    currentV = add_UF(currentV, mul_UF(convertInteger_UF(1914210), power_UF(arg, 8)));
    currentV = add_UF(currentV, mul_UF(convertInteger_UF(4744640), power_UF(arg, 6)));
    currentV = add_UF(currentV, mul_UF(convertInteger_UF(1891200), power_UF(arg, 4)));
    currentV = add_UF(currentV, mul_UF(convertInteger_UF(78720), power_UF(arg, 2)));
    currentV = add_UF(currentV, convertInteger_UF(256));
    currentV = div_UF(currentV, mul_UF(convertInteger_UF(322560), power_UF(convertInteger_UF(n), 5)));
    currentV = div_UF(currentV, power_UF(aux2, 7));
    currentV = div_UF(currentV, aux1);
    currentV = sub_UF(currentV, div_UF(convertInteger_UF(1), mul_UF(convertInteger_UF(1260), power_UF(convertInteger_UF(n), 5))));
    sumV = add_UF(sumV, currentV);

    //V6
    currentV = convertInteger_UF(48);
    currentV = add_UF(currentV, mul_UF(convertInteger_UF(2580), power_UF(arg, 2)));
    currentV = add_UF(currentV, mul_UF(convertInteger_UF(14884), power_UF(arg, 4)));
    currentV = add_UF(currentV, mul_UF(convertInteger_UF(17493), power_UF(arg, 6)));
    currentV = add_UF(currentV, mul_UF(convertInteger_UF(4242), power_UF(arg, 8)));
    currentV = add_UF(currentV, mul_UF(convertInteger_UF(103), power_UF(arg, 10)));
    currentV = mul_UF(currentV, power_UF(arg, 2));
    currentV = div_UF(currentV, mul_UF(convertInteger_UF(192), power_UF(convertInteger_UF(n), 6)));
    currentV = div_UF(currentV, power_UF(aux2, 9));
    sumV = add_UF(sumV, currentV);

    //V7
    currentV = mul_UF(convertInteger_UF(881664), power_UF(arg, 2));
    currentV = add_UF(currentV, mul_UF(convertInteger_UF(99783936), power_UF(arg, 4)));
    currentV = add_UF(currentV, mul_UF(convertInteger_UF(1135145088), power_UF(arg, 6)));
    currentV = add_UF(currentV, mul_UF(convertInteger_UF(2000000000), power_UF(arg, 8)));
    currentV = add_UF(currentV, mul_UF(convertInteger_UF(884531440), power_UF(arg, 8)));
    currentV = add_UF(currentV, mul_UF(convertInteger_UF(1965889800), power_UF(arg, 10)));
    currentV = add_UF(currentV, mul_UF(convertInteger_UF(318291750), power_UF(arg, 12)));
    currentV = add_UF(currentV, mul_UF(convertInteger_UF(5635995), power_UF(arg, 14)));
    currentV = sub_UF(currentV, convertInteger_UF(2048));
    currentV = div_UF(currentV, mul_UF(convertInteger_UF(3440640), power_UF(convertInteger_UF(n), 7)));
    currentV = div_UF(currentV, power_UF(aux2, 10));
    currentV = div_UF(currentV, aux1);
    currentV = add_UF(currentV, div_UF(convertInteger_UF(1), mul_UF(convertInteger_UF(1680), power_UF(convertInteger_UF(n), 7))));
    sumV = add_UF(sumV, currentV);

    //V8
    currentV = convertInteger_UF(1024);
    currentV = add_UF(currentV, mul_UF(convertInteger_UF(248320), power_UF(arg, 2)));
    currentV = add_UF(currentV, mul_UF(convertInteger_UF(5095936), power_UF(arg, 4)));
    currentV = add_UF(currentV, mul_UF(convertInteger_UF(24059968), power_UF(arg, 6)));
    currentV = add_UF(currentV, mul_UF(convertInteger_UF(34280896), power_UF(arg, 8)));
    currentV = add_UF(currentV, mul_UF(convertInteger_UF(15252048), power_UF(arg, 10)));
    currentV = add_UF(currentV, mul_UF(convertInteger_UF(1765936), power_UF(arg, 12)));
    currentV = add_UF(currentV, mul_UF(convertInteger_UF(23797), power_UF(arg, 14)));
    currentV = mul_UF(currentV, power_UF(arg, 2));
    currentV = div_UF(currentV, mul_UF(convertInteger_UF(4096), power_UF(convertInteger_UF(n), 8)));
    currentV = div_UF(currentV, power_UF(aux2, 12));
    sumV = add_UF(sumV, currentV);

    //////compute result
    result = div_UF(convertInteger_UF(1), sqrt_UF(aux1));
    current1 = log_UF(div_UF(arg, add_UF(convertInteger_UF(1), aux1)));
    current1 = add_UF(current1, aux1);
    current1 = mul_UF(current1, convertInteger_UF(n));

    current2 = add_UF(current2, sumV);
    current2 = add_UF(current2, factorialSeries_UF(n));

    result = mul_UF(result, exp_UF(sub_UF(current2, current1)));
    result = div_UF(result, sqrt_UF(Pi));
    result = mul_UF(result, sqrt_UF(convertInteger_UF(2)));
    result = div_UF(result, sqrt_UF(convertInteger_UF(n)));
    result->sign = -1;
    return result;
}

//////////////////////////////////////////////////////////////////////
//                        jnMeisselSecond_UF                        //
//////////////////////////////////////////////////////////////////////
Unifloat* jnMeisselSecond_UF(int n, Unifloat* x)
{
    Unifloat* result, * current;
    Unifloat* P = createZero_UF();
    Unifloat* Q = createZero_UF();
    Unifloat* arg, * cot, * sec;

    arg = acos_UF(div_UF(convertInteger_UF(n), x));
    cot = div_UF(convertInteger_UF(1), tan_UF(arg));
    sec = div_UF(x, convertInteger_UF(n));

    //P1
    current = mul_UF(convertInteger_UF(4), power_UF(sec, 2));
    current = add_UF(current, power_UF(sec, 4));
    current = mul_UF(current, power_UF(cot, 6));
    current = div_UF(current, mul_UF(convertInteger_UF(16), power_UF(convertInteger_UF(n), 2)));
    P = add_UF(P, current);

    //P2
    current = mul_UF(convertInteger_UF(32), power_UF(sec, 2));
    current = add_UF(current, mul_UF(convertInteger_UF(288), power_UF(sec, 4)));
    current = add_UF(current, mul_UF(convertInteger_UF(232), power_UF(sec, 6)));
    current = add_UF(current, mul_UF(convertInteger_UF(13), power_UF(sec, 8)));
    current = mul_UF(current, power_UF(cot, 12));
    current = div_UF(current, mul_UF(convertInteger_UF(128), power_UF(convertInteger_UF(n), 4)));
    P = sub_UF(P, current);

    //P3
    current = mul_UF(convertInteger_UF(48), power_UF(sec, 2));
    current = add_UF(current, mul_UF(convertInteger_UF(2580), power_UF(sec, 4)));
    current = add_UF(current, mul_UF(convertInteger_UF(14884), power_UF(sec, 6)));
    current = add_UF(current, mul_UF(convertInteger_UF(17493), power_UF(sec, 8)));
    current = add_UF(current, mul_UF(convertInteger_UF(4242), power_UF(sec, 10)));
    current = add_UF(current, mul_UF(convertInteger_UF(103), power_UF(sec, 12)));
    current = mul_UF(current, power_UF(cot, 18));
    current = div_UF(current, mul_UF(convertInteger_UF(192), power_UF(convertInteger_UF(n), 6)));
    P = add_UF(P, current);

    //P4
    current = convertInteger_UF(1024);
    current = add_UF(current, mul_UF(convertInteger_UF(248320), power_UF(sec, 2)));
    current = add_UF(current, mul_UF(convertInteger_UF(5095936), power_UF(sec, 4)));
    current = add_UF(current, mul_UF(convertInteger_UF(24059968), power_UF(sec, 6)));
    current = add_UF(current, mul_UF(convertInteger_UF(34280896), power_UF(sec, 8)));
    current = add_UF(current, mul_UF(convertInteger_UF(15252048), power_UF(sec, 10)));
    current = add_UF(current, mul_UF(convertInteger_UF(1765936), power_UF(sec, 12)));
    current = add_UF(current, mul_UF(convertInteger_UF(23797), power_UF(sec, 14)));
    current = mul_UF(current, power_UF(cot, 24));
    current = mul_UF(current, power_UF(sec, 2));
    current = div_UF(current, mul_UF(convertInteger_UF(4096), power_UF(convertInteger_UF(n), 8)));
    P = sub_UF(P, current);

    //Q1
    current = add_UF(convertInteger_UF(2), mul_UF(convertInteger_UF(3), power_UF(sec, 2)));
    current = mul_UF(current, power_UF(cot, 3));
    current = div_UF(current, mul_UF(convertInteger_UF(24), convertInteger_UF(n)));
    current->sign = (-1)*current->sign;
    current = add_UF(current, mul_UF(convertInteger_UF(n), sub_UF(tan_UF(arg), arg)));
    Q = add_UF(Q, current);

    //Q2
    current = convertInteger_UF(16);
    current = sub_UF(current, mul_UF(convertInteger_UF(1512), power_UF(sec, 2)));
    current = sub_UF(current, mul_UF(convertInteger_UF(3654), power_UF(sec, 4)));
    current = sub_UF(current, mul_UF(convertInteger_UF(375), power_UF(sec, 6)));
    current = mul_UF(current, power_UF(cot, 9));
    current = div_UF(current, mul_UF(convertInteger_UF(5760), power_UF(convertInteger_UF(n), 3)));
    Q = sub_UF(Q, current);

    //Q3
    current = convertInteger_UF(256);
    current = add_UF(current, mul_UF(convertInteger_UF(78720), power_UF(sec, 2)));
    current = add_UF(current, mul_UF(convertInteger_UF(1891200), power_UF(sec, 4)));
    current = add_UF(current, mul_UF(convertInteger_UF(4744640), power_UF(sec, 6)));
    current = add_UF(current, mul_UF(convertInteger_UF(1914210), power_UF(sec, 8)));
    current = add_UF(current, mul_UF(convertInteger_UF(67599), power_UF(sec, 10)));
    current = mul_UF(current, power_UF(cot, 15));
    current = div_UF(current, mul_UF(convertInteger_UF(322560), power_UF(convertInteger_UF(n), 5)));
    Q = sub_UF(Q, current);

    //Q4
    current = mul_UF(convertInteger_UF(881664), power_UF(sec, 2));
    current = add_UF(current, mul_UF(convertInteger_UF(99783936), power_UF(sec, 4)));
    current = add_UF(current, mul_UF(convertInteger_UF(1135145088), power_UF(sec, 6)));
    current = add_UF(current, mul_UF(convertInteger_UF(2000000000), power_UF(sec, 8)));
    current = add_UF(current, mul_UF(convertInteger_UF(884531440), power_UF(sec, 8)));
    current = add_UF(current, mul_UF(convertInteger_UF(1965889800), power_UF(sec, 10)));
    current = add_UF(current, mul_UF(convertInteger_UF(318291750), power_UF(sec, 12)));
    current = add_UF(current, mul_UF(convertInteger_UF(5635995), power_UF(sec, 14)));
    current = sub_UF(current, convertInteger_UF(2048));
    current = mul_UF(current, power_UF(cot, 21));
    current = div_UF(current, mul_UF(convertInteger_UF(3440640), power_UF(convertInteger_UF(n), 7)));
    Q = add_UF(Q, current);

    //////result//////
    result = sqrt_UF(div_UF(mul_UF(convertInteger_UF(2), cot), mul_UF(convertInteger_UF(n), Pi)));
    result = div_UF(result, exp_UF(P));
    
    result = mul_UF(result, cos_UF(sub_UF(Q, div_UF(Pi, convertInteger_UF(4)))));

    return result;
}

//////////////////////////////////////////////////////////////////////
//                        ynMeisselSecond_UF                        //
//////////////////////////////////////////////////////////////////////
Unifloat* ynMeisselSecond_UF(int n, Unifloat* x)
{
    Unifloat* result, * current;
    Unifloat* P = createZero_UF();
    Unifloat* Q = createZero_UF();
    Unifloat* arg, * cot, * sec;

    arg = acos_UF(div_UF(convertInteger_UF(n), x));
    cot = div_UF(convertInteger_UF(1), tan_UF(arg));
    sec = div_UF(x, convertInteger_UF(n));

    //P1
    current = mul_UF(convertInteger_UF(4), power_UF(sec, 2));
    current = add_UF(current, power_UF(sec, 4));
    current = mul_UF(current, power_UF(cot, 6));
    current = div_UF(current, mul_UF(convertInteger_UF(16), power_UF(convertInteger_UF(n), 2)));
    P = add_UF(P, current);

    //P2
    current = mul_UF(convertInteger_UF(32), power_UF(sec, 2));
    current = add_UF(current, mul_UF(convertInteger_UF(288), power_UF(sec, 4)));
    current = add_UF(current, mul_UF(convertInteger_UF(232), power_UF(sec, 6)));
    current = add_UF(current, mul_UF(convertInteger_UF(13), power_UF(sec, 8)));
    current = mul_UF(current, power_UF(cot, 12));
    current = div_UF(current, mul_UF(convertInteger_UF(128), power_UF(convertInteger_UF(n), 4)));
    P = sub_UF(P, current);

    //P3
    current = mul_UF(convertInteger_UF(48), power_UF(sec, 2));
    current = add_UF(current, mul_UF(convertInteger_UF(2580), power_UF(sec, 4)));
    current = add_UF(current, mul_UF(convertInteger_UF(14884), power_UF(sec, 6)));
    current = add_UF(current, mul_UF(convertInteger_UF(17493), power_UF(sec, 8)));
    current = add_UF(current, mul_UF(convertInteger_UF(4242), power_UF(sec, 10)));
    current = add_UF(current, mul_UF(convertInteger_UF(103), power_UF(sec, 12)));
    current = mul_UF(current, power_UF(cot, 18));
    current = div_UF(current, mul_UF(convertInteger_UF(192), power_UF(convertInteger_UF(n), 6)));
    P = add_UF(P, current);

    //P4
    current = convertInteger_UF(1024);
    current = add_UF(current, mul_UF(convertInteger_UF(248320), power_UF(sec, 2)));
    current = add_UF(current, mul_UF(convertInteger_UF(5095936), power_UF(sec, 4)));
    current = add_UF(current, mul_UF(convertInteger_UF(24059968), power_UF(sec, 6)));
    current = add_UF(current, mul_UF(convertInteger_UF(34280896), power_UF(sec, 8)));
    current = add_UF(current, mul_UF(convertInteger_UF(15252048), power_UF(sec, 10)));
    current = add_UF(current, mul_UF(convertInteger_UF(1765936), power_UF(sec, 12)));
    current = add_UF(current, mul_UF(convertInteger_UF(23797), power_UF(sec, 14)));
    current = mul_UF(current, power_UF(cot, 24));
    current = mul_UF(current, power_UF(sec, 2));
    current = div_UF(current, mul_UF(convertInteger_UF(4096), power_UF(convertInteger_UF(n), 8)));
    P = sub_UF(P, current);

    //Q1
    current = add_UF(convertInteger_UF(2), mul_UF(convertInteger_UF(3), power_UF(sec, 2)));
    current = mul_UF(current, power_UF(cot, 3));
    current = div_UF(current, mul_UF(convertInteger_UF(24), convertInteger_UF(n)));
    current->sign = (-1)*current->sign;
    current = add_UF(current, mul_UF(convertInteger_UF(n), sub_UF(tan_UF(arg), arg)));
    Q = add_UF(Q, current);

    //Q2
    current = convertInteger_UF(16);
    current = sub_UF(current, mul_UF(convertInteger_UF(1512), power_UF(sec, 2)));
    current = sub_UF(current, mul_UF(convertInteger_UF(3654), power_UF(sec, 4)));
    current = sub_UF(current, mul_UF(convertInteger_UF(375), power_UF(sec, 6)));
    current = mul_UF(current, power_UF(cot, 9));
    current = div_UF(current, mul_UF(convertInteger_UF(5760), power_UF(convertInteger_UF(n), 3)));
    Q = sub_UF(Q, current);

    //Q3
    current = convertInteger_UF(256);
    current = add_UF(current, mul_UF(convertInteger_UF(78720), power_UF(sec, 2)));
    current = add_UF(current, mul_UF(convertInteger_UF(1891200), power_UF(sec, 4)));
    current = add_UF(current, mul_UF(convertInteger_UF(4744640), power_UF(sec, 6)));
    current = add_UF(current, mul_UF(convertInteger_UF(1914210), power_UF(sec, 8)));
    current = add_UF(current, mul_UF(convertInteger_UF(67599), power_UF(sec, 10)));
    current = mul_UF(current, power_UF(cot, 15));
    current = div_UF(current, mul_UF(convertInteger_UF(322560), power_UF(convertInteger_UF(n), 5)));
    Q = sub_UF(Q, current);

    //Q4
    current = mul_UF(convertInteger_UF(881664), power_UF(sec, 2));
    current = add_UF(current, mul_UF(convertInteger_UF(99783936), power_UF(sec, 4)));
    current = add_UF(current, mul_UF(convertInteger_UF(1135145088), power_UF(sec, 6)));
    current = add_UF(current, mul_UF(convertInteger_UF(2000000000), power_UF(sec, 8)));
    current = add_UF(current, mul_UF(convertInteger_UF(884531440), power_UF(sec, 8)));
    current = add_UF(current, mul_UF(convertInteger_UF(1965889800), power_UF(sec, 10)));
    current = add_UF(current, mul_UF(convertInteger_UF(318291750), power_UF(sec, 12)));
    current = add_UF(current, mul_UF(convertInteger_UF(5635995), power_UF(sec, 14)));
    current = sub_UF(current, convertInteger_UF(2048));
    current = mul_UF(current, power_UF(cot, 21));
    current = div_UF(current, mul_UF(convertInteger_UF(3440640), power_UF(convertInteger_UF(n), 7)));
    Q = add_UF(Q, current);

    //////result//////
    result = sqrt_UF(div_UF(mul_UF(convertInteger_UF(2), cot), mul_UF(convertInteger_UF(n), Pi)));
    result = div_UF(result, exp_UF(P));
    result = mul_UF(result, sin_UF(sub_UF(Q, div_UF(Pi, convertInteger_UF(4)))));
    return result;
}


//////////////////////////////////////////////////////////////////////
//                         jnRecurrent_UF                           //
//////////////////////////////////////////////////////////////////////
Unifloat* jnRecurrent_UF(int n, Unifloat* x)
{
    int i, NUM = 1;
    Unifloat* INDEX = convertInteger_UF(n);
    Unifloat* aux, * tmp, * A1, * B1, * A2, * B2, * CF, * J0, * Q0, * Q1;
    Unifloat* One, * Two, * Zero;
    Zero = createZero_UF();
    One = createOne_UF();
    Two = clone(One);
    Two->exp++;

    J0 = j0_UF(x);
    A1 = clone(J0);
    B1 = j1_UF(x);
    if(compare_UF(abs_UF(x), abs_UF(INDEX)) != -1)
    {
        for(i = 1; i < n; i++)
        {
            aux = clone(B1);
            B1 = mul_UF(B1, convertInteger_UF(i));
            B1 = div_UF(B1, x);
            B1->exp++;
            B1 = sub_UF(B1, A1);
            A1 = clone(aux);
        }
        return B1;
    }
    else
    {
        //////determine NUM///////
        B1 = div_UF(createOne_UF(), x);
        B1->exp++;
        A1 = mul_UF(INDEX, B1);
        Q0 = clone(A1);
        Q1 = mul_UF(Q0, add_UF(A1, B1));
        Q1 = sub_UF(Q1, One);
        aux = add_UF(A1, B1);
        while(Q1->exp < 50)
        {
            NUM += 1;
            aux = add_UF(aux, B1);
            tmp = mul_UF(aux, Q1);
            tmp = sub_UF(tmp, Q0);
            Q0 = clone(Q1);
            Q1 = clone(tmp);
        }

        //////compute CF//////
        A1 = clone(One);
        A2 = clone(Zero);
        
        B1 = clone(Zero);
        B2 = clone(One);

        for(i = 0; i < NUM; i++)
        {
            aux = clone(A2);
            A2 = sub_UF(div_UF(mul_UF(mul_UF(Two, A2), 
               add_UF(INDEX, convertInteger_UF(i))), x), A1);
            A1 = clone(aux);
            aux = clone(B2);
            B2 = sub_UF(div_UF(mul_UF(mul_UF(Two, B2), 
               add_UF(INDEX, convertInteger_UF(i))), x), B1);
            B1 = clone(aux);
        }
        
        CF = div_UF(A2, B2);
        CF->sign = (-1) * CF->sign;

        //////backward recurrence//////
        A1 = clone(CF);
        B1 = createOne_UF();
        for(i = n - 1; i > 0; i--)
        {
            aux = clone(B1);
            B1 = mul_UF(B1, convertInteger_UF(i));
            B1 = div_UF(B1, x);
            B1->exp++;
            B1 = sub_UF(B1, A1);
            A1 = clone(aux);
        }
        B1 = div_UF(mul_UF(CF, J0), B1);
        return B1;
    }
}

//////////////////////////////////////////////////////////////////////
//                         ynRecurrent_UF                           //
//////////////////////////////////////////////////////////////////////
Unifloat* ynRecurrent_UF(int n, Unifloat* x)
{
    int i;
    Unifloat* aux, * A1, * B1;

    A1 = y0_UF(x);
    B1 = y1_UF(x);
    for(i = 1; i < n; i++)
    {
        aux = clone(B1);
        B1 = mul_UF(B1, convertInteger_UF(i));
        B1 = div_UF(B1, x);
        B1->exp++;
        B1 = sub_UF(B1, A1);
        A1 = clone(aux);
    }
    return B1;
}
