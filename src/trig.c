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
#include "unifloat/trig.h"
#include "unifloat/exp.h"
#include "unifloat/debug.h"

//////////////////////////////////////////////////////////////////////////
//                          sin_UF                                      //
//////////////////////////////////////////////////////////////////////////
Unifloat* sin_UF(Unifloat* x)
{
    /*
    res - contains result
    type - contains current type
    pi, One, Two, Four - are constant numbers
    arg - is the same as x, using for x don't change in calculating
    tmp - auxiliary variable
    eps - contains EPSILON
    */
    Unifloat* res, *optmp;
    Unifloat* One;
    Unifloat* Two;
    Unifloat* Four;
    Unifloat* arg;
    Unifloat* tmp;

    arg = clone(x);

    if(isNan_UF(arg)) return arg; // atan(NaN)=NaN without error

    if(isInfinity_UF(arg))
    {    
        // sin(Infinity) = NaN
        arg->kind = NaN;
        return arg;
    }

    if (arg->sign == -1)
    {   
        // sin(-x) = -sin(x);
        arg->sign = 1;
        res = sin_UF(arg);
        res->sign = (-1)*res->sign;
        delete_UF(arg);
        return res;
    }

    arg = call1_arg1(&arrangeArgument_UF, arg);

    if (compare_UF(arg, Pi) == 1)
    {
        optmp = sub_UF(arg, Pi);
        res = sin_UF(optmp);
        delete_UF(optmp);
        res->sign = (-1)*res->sign;
        delete_UF(arg);
        return res;
    }

    Two = convertInteger_UF(2);
    optmp = div_UF(Pi, Two);
    if(compare_UF(arg, optmp) == 1) 
    {
        arg = call2_arg2(&sub_UF, Pi, arg);                    
        if (arg->exp <= -PRECISION) 
        {
            delete_UF(arg);
            arg = createZero_UF();
        }
    }
    delete_UF(optmp);

    Four = convertInteger_UF(4);
    tmp = div_UF(Pi, Four);
    delete_UF(Four);

    if(compareWithPrecision_UF(arg, tmp, PRECISION-2) == 1)
    {
        optmp = div_UF(Pi, Two);
        arg = call2_arg2(&sub_UF, optmp, arg);
    
        res = cos_UF(arg);
        delete_UFs(arg, Two, tmp, optmp, NULL);
        return res;
    }

    if(isZero_UF(arg))
    {
        delete_UFs(arg, Two, tmp, NULL);
        return createZero_UF();
    }

    One = convertInteger_UF(1);

    optmp = div_UF(arg, Two);
    res = tan_UF(optmp);
    delete_UF(optmp);

    optmp = mul_UF(res, res);
    delete_UF(tmp);
    tmp = add_UF(One, optmp);
    delete_UF(optmp);

    optmp = mul_UF(Two, res);
    delete_UF(res);
    res = div_UF(optmp, tmp);
    delete_UFs(optmp, tmp, NULL);

    optmp = abs_UF(res);
    if(compare_UF(optmp, One) == 1)
    {
        if(res->sign == 1)   
        {
            delete_UF(res);         
            res = createOne_UF();
        }
        else 
        {
            delete_UF(res);
            res = createOne_UF();
            res->sign =  - 1;
        }
    }
    delete_UFs(optmp, One, Two, arg, NULL);
    return res;
}

//////////////////////////////////////////////////////////////////////////
//                          cos_UF                                      //
//////////////////////////////////////////////////////////////////////////
Unifloat* cos_UF(Unifloat* x)
{
    /*
    res - contains result
    pi, One, Two, Four - are constant numbers
    arg - is the same as x, using for x don't change in calculating
    tmp - auxiliary variable
    eps - contains EPSILON
    */
    Unifloat* res, *optmp, *optmp2, *optmp3;
    //Unifloat* eps = clone(Epsilon);
    Unifloat* One;
    Unifloat* Two;
    Unifloat* Four;
    Unifloat* arg;
    Unifloat* tmp;

    arg = clone(x);

    if(isNan_UF(arg)) return arg; // atan(NaN)=NaN without 
    
    if(isInfinity_UF(arg))
    { // cos(Infinity) = NaN
        arg->kind = NaN;
        return arg;
    }

    if (arg->sign == -1)
    {
        arg->sign = 1;
        optmp = cos_UF(arg);
        delete_UF(arg);
        return optmp;
    }

    arg = call1_arg1(&arrangeArgument_UF, arg);

    if(compare_UF(arg, Pi) == 1)
    { // if arg>pi then res=-cos(arg-pi)
        optmp = sub_UF(arg, Pi);
        res = cos_UF(optmp);
        res->sign = (-1)*res->sign;
        delete_UFs(optmp, arg, NULL);
        return res;
    }

    Two = convertInteger_UF(2);
    optmp2 = div_UF(Pi, Two);
    if(compare_UF(arg, optmp2) == 1) 
    { // if arg>pi/2 then res=-cos(pi-arg)
        optmp = sub_UF(Pi, arg);
        res = cos_UF(optmp);
        res->sign = (-1)*res->sign;
        delete_UFs(optmp, arg, Two, optmp2, NULL);
        return res;
    }
    delete_UF(optmp2);
    Four = convertInteger_UF(4);
    tmp = div_UF(Pi, Four);
    delete_UF(Four);
    if(compareWithPrecision_UF(arg, tmp, PRECISION-2) == 1)
    {
        optmp = div_UF(Pi, Two);
        optmp2 = sub_UF(optmp, arg);
        optmp3 = sin_UF(optmp2);
        delete_UFs(optmp, optmp2, tmp, arg, Two, NULL);
        return optmp3;
    }

//    if(compare_UF(abs_UF(arg), eps) != 1) // cos(arg)=1, |arg|<eps
//        return createOne_UF();

    optmp = div_UF(arg, Two);
    res = tan_UF(optmp);
    delete_UF(optmp);

    delete_UF(tmp);
    tmp = mul_UF(res, res);

    One = createOne_UF();
    optmp = sub_UF(One, tmp);
    optmp2 = add_UF(One, tmp);
    delete_UF(res);
    res = div_UF(optmp, optmp2);
    delete_UFs(optmp, optmp2, NULL);

    optmp = abs_UF(res);
    if(compare_UF(optmp, One) == 1)
    {
        if(res->sign == 1) // if res>1 then return 1
            res = createOne_UF();
        else 
        { // if res<-1 then return -1
            res = createOne_UF();
            res->sign =  - 1;
        }
    }
    delete_UFs(One, Two, arg, tmp, optmp, NULL);
    return res;
}

//////////////////////////////////////////////////////////////////////////
//                          tan_UF                                      //
//////////////////////////////////////////////////////////////////////////
Unifloat* tan_UF(Unifloat* x)
{
    /*
    res - contains result
    pi, One, Two, Four - are constant numbers
    arg - is the same as x, using for x don't change in calculating
    tmp - auxiliary variable
    eps - contains EPSILON
    */
    Bool cmpres;
    Unifloat* res, *optmp, *optmp2;
    Unifloat* eps;
    Unifloat* One;
    Unifloat* Two;
    Unifloat* Four;
    Unifloat* arg, * tmp;

    arg = clone(x);
    
    if(isNan_UF(arg)) return arg; // tan(NaN)=NaN without error

    if(isInfinity_UF(arg))
    { // tan(Infinity) = NaN
        arg->kind = NaN;
        return arg;
    }

    if (arg->sign == -1)
    {
        arg->sign = 1;
        res = tan_UF(arg);
        res->sign = -res->sign;
        delete_UF(arg);
        return res;
    }

    arg = call1_arg1(&arrangeArgument_UF, arg);

    if(compare_UF(arg, Pi) == 1) // if arg>pi then arg:=arg-pi
        arg = call2_arg1(&sub_UF, arg, Pi);

/*  tan(arg)=inf, |arg-pi/2|<eps
    tmp = abs_UF(sub_UF(arg, div_UF(pi, Two)));
    if(compare_UF(tmp, eps) != 1)
    {    
        arg = createZero_UF(BASE, PRECISION);
        res->kind = Infinity;
        return res;
    }
*/

    Two = convertInteger_UF(2);
    optmp = div_UF(Pi, Two);
    cmpres = compare_UF(arg,optmp );
    delete_UF(optmp);

    if(cmpres == 1)
    { // if arg>pi/2 then res=-tan(pi-arg)
        optmp = sub_UF(Pi, arg);
        res = tan_UF(optmp);
        res->sign = (-1)*res->sign;
        delete_UFs(optmp, arg, Two, NULL);
        return res;
    }

    One = createOne_UF();
    Four = convertInteger_UF(4);

    optmp = div_UF(Pi, Four);
    cmpres = compare_UF(arg, optmp);
    delete_UF(optmp);

    if(cmpres == 1)
    { // if arg>pi/4 then res=1/tan(pi/2-arg)
        optmp = div_UF(Pi, Two);
        optmp2 = sub_UF(optmp, arg);
        res = tan_UF(optmp2);
        res = call2_arg2(&div_UF, One, res);
        delete_UFs(optmp2, optmp, arg, One, Two, Four, NULL);
        return res;
    }

    if(isZero_UF(arg)) // tan(0)=0
    {
        delete_UFs(arg, One, Two, Four, NULL);
        return createZero_UF();
    }

    eps = clone(Epsilon);
    optmp = abs_UF(arg);
    if(compare_UF(optmp, eps) != 1) // tan(arg)=arg, |arg|<eps
    {
        delete_UFs(One, Two, Four, eps, optmp, NULL);
        return arg;
    }
    delete_UF(optmp);

    optmp = div_UF(arg, Two);
    res = calcTan_UF(optmp, 1);
    delete_UF(optmp);

    optmp = mul_UF(res, res);
    tmp = sub_UF(One, optmp);
    delete_UF(optmp);

    optmp = mul_UF(Two, res);
    delete_UF(res);
    res = div_UF(optmp, tmp);
    delete_UFs(optmp, tmp, arg, One, Two, Four, eps, NULL);
    return res;
}

//////////////////////////////////////////////////////////////////////////
//                          asin_UF                                     //
//////////////////////////////////////////////////////////////////////////
Unifloat* asin_UF(Unifloat* x)
{
    /*
    res - contains result
    pi, One, Two, Four - are constant numbers
    arg - is the same as x, using for x don't change in calculating
    tmp - auxiliary variable
    eps - contains EPSILON
    */
    Unifloat* res, * tmp, *optmp, *optmp2;
    Unifloat* eps = createOne_UF();
    Unifloat* One = createOne_UF();
    Unifloat* Two;
    Unifloat* arg;
    
    eps->exp = -digMant_DoubleT;

    arg = clone(x);

    if(isNan_UF(arg)) 
    {
        delete_UFs(One, eps, NULL);
        return arg; // asin(NaN)=NaN without error
    }

    optmp = abs_UF(arg);
    if((compare_UF(optmp, One) == 1)
        ||(isInfinity_UF(arg)))
    { // asin(agr<-1||arg>1||NaN||Infinity)==NaN
        arg->kind = NaN;
        delete_UFs(One, eps, optmp, NULL);
        return arg;
    }
    delete_UF(optmp);

    if(isZero_UF(arg)) // asin(0)=0
    {
        delete_UFs(One, eps, NULL);
        return arg;
    }

    optmp = abs_UF(arg);
    if(compare_UF(optmp, eps) != 1) //asin(arg)=arg, |arg|<eps
    {
        delete_UFs(One, eps, optmp, NULL);
        return arg;
    }
    delete_UF(optmp);

    if (arg->sign == -1)
    {
        arg->sign = 1;
        res = asin_UF(arg);
        res->sign = -res->sign;
        delete_UFs(arg, eps, One, NULL);
        return res;
    }

    Two = convertInteger_UF(2);
    optmp = sub_UF(arg, One);
    optmp = call1_arg1(&abs_UF,optmp);
    if(compare_UF(optmp, eps) != 1)
    {
        delete_UF(optmp);
        optmp = div_UF(Pi, Two);
        delete_UFs(arg, Two, One, eps, NULL);
        return optmp; // asin(arg)=pi/2, |arg-1|<eps
    }
    delete_UF(optmp);

    optmp = add_UF(arg, One);
    optmp = call1_arg1(&abs_UF, optmp);

    if(compare_UF(optmp, eps) != 1)
    { // asin(arg)=-pi/2, |arg+1|<eps
        res = div_UF(Pi, Two);
        res->sign = (-1)*res->sign;
        delete_UFs(optmp, arg, Two, One, eps);
        return res;
    }
    delete_UF(optmp);

    optmp = mul_UF(arg, arg);
    optmp2 = sub_UF(One, optmp);
    res = sqrt_UF(optmp2);
    delete_UFs(optmp, optmp2, NULL);

    optmp = div_UF(arg, res);
    delete_UF(res);
    res = atan_UF(optmp);
    delete_UF(optmp);

    optmp = div_UF(Pi, Two);
    if(compare_UF(res, optmp) == 1)
    {
        delete_UF(res);
        res = div_UF(Pi, Two); // if res>pi/2 then return pi/2
    }
    delete_UF(optmp);

    tmp = div_UF(Pi, Two);
    tmp->sign = (-1)*tmp->sign;
    if(compare_UF(tmp, res) == 1)
    { // if res<-pi/2 then return -pi/2
        delete_UF(res);
        res = div_UF(Pi, Two);
        res->sign = (-1)*res->sign;
    }

    delete_UFs(eps, One, Two, arg, tmp, NULL);
    return res;
}

//////////////////////////////////////////////////////////////////////////
//                          acos_UF                                     //
//////////////////////////////////////////////////////////////////////////
Unifloat* acos_UF(Unifloat* x)
{
    /*
    res - contains result
    pi, One, Two, Zero - are constant numbers
    arg - is the same as x, using for x don't change in calculating
    eps - contains EPSILON
    */
    Unifloat* res, *optmp, *optmp2;
    Unifloat* Zero = createZero_UF();
    Unifloat* One = createOne_UF();
    Unifloat* Two;
    Unifloat* arg;

    arg = clone(x);

    if(isNan_UF(arg)) 
    {
        delete_UFs(Zero, One, NULL);
        return arg; // acos(NaN)=NaN without error
    }

    optmp = abs_UF(arg);
    if((compare_UF(optmp, One) == 1)//
        ||(isInfinity_UF(arg)))
    { // acos(agr<-1||arg>1||Infinity)==NaN
        arg->kind = NaN;
        delete_UFs(optmp, One, Zero, NULL);
        return arg;
    }
    delete_UF(optmp);

    if(compare_UF(arg, One) == 0) // acos(1)==0
    {
        delete_UFs(arg, One, NULL);
        return Zero;
    }

    Two = convertInteger_UF(2);
    if(isZero_UF(arg)) // acos(0)=pi/2
    {
        optmp = div_UF(Pi, Two);
        delete_UFs(arg, Two, One, Zero, NULL);
        return optmp;
    }

    optmp = div_UF(Pi, Two);
    optmp2 = asin_UF(arg);
    res = sub_UF(optmp, optmp2);
    delete_UFs(optmp, optmp2, NULL);
    if(compare_UF(res, Pi) == 1) 
    {
        delete_UF(res);
        res = Pi; // if res>pi then return pi
    }
    
    if(compare_UF(Zero, res) == 1)
        res = Zero; // if res<0 then return 0

    delete_UFs(arg, Two, One, Zero, NULL);
    return res;
}

//////////////////////////////////////////////////////////////////////////
//                          atan_UF                                     //
//////////////////////////////////////////////////////////////////////////
Unifloat* atan_UF(Unifloat* x)
{
    /*
    res - contains result
    pi, One, Two, SqrtThree, Six - are constant numbers
    arg - is the same as x, using for x don't change in calculating
    tmp - auxiliary variable
    eps - contains EPSILON
    new_arg - auxiliary variable for calcalating
    */
    Unifloat *optmp, *optmp2;
    Unifloat* eps;
    Unifloat* One;
    Unifloat* Two;
    Unifloat* SqrtThree;
    Unifloat* Six;
    Unifloat* new_arg, * res, *tmp, *arg;

    arg = clone(x);

    if(isNan_UF(arg)) 
        return arg; // atan(NaN)=NaN without error

    Two = convertInteger_UF(2);
    if(isInfinity_UF(arg)) 
    { // atan(inf)=pi/2;
        res=div_UF(Pi, Two);
        res->sign=arg->sign;
        delete_UFs(arg, Two, NULL);
        return res;
    }

    if(arg->sign ==  -1)
    { // if arg<0 then res=-atan(-arg)
        arg->sign=(-1)*arg->sign;
        res = atan_UF(arg);
        res->sign = (-1)*res->sign;
        delete_UFs(arg, Two, NULL);
        return res;
    }

    One = createOne_UF();
    if(compare_UF(arg, One) == 1)
    { // if arg>1 then res=pi/2-atan(1/arg)
        optmp = div_UF(One, arg);
        tmp = atan_UF(optmp);
        delete_UF(optmp);

        optmp = div_UF(Pi, Two);
        optmp = call2_arg1(&sub_UF, optmp, tmp);
        delete_UFs(tmp, arg, One, Two, NULL);
        return optmp;
    }

    if(isZero_UF(arg)) // atan(0)==0
    {
        delete_UFs(arg, One, Two, NULL);
        return createZero_UF();
    }

    eps = createOne_UF();
    eps->exp = -digMant_DoubleT;
    optmp = abs_UF(arg);
    if(compare_UF(optmp, eps) != 1) // atan(arg)=arg, |arg|<eps
    {
        delete_UFs(optmp, arg, eps, One, Two, NULL);
        return arg;
    }
    delete_UF(optmp);

    optmp = convertInteger_UF(3);
    SqrtThree = sqrt_UF(optmp);
    delete_UF(optmp);
    optmp2 = sub_UF(Two, SqrtThree);
    if(compare_UF(arg, optmp2) == 1) //x > 2- sqrt(3)
    {
        optmp = mul_UF(arg, SqrtThree);
        tmp = sub_UF(optmp, One); // x*sqrt(3) - 1
        delete_UF(optmp);
        
        optmp = add_UF(arg, SqrtThree);
        new_arg = div_UF(tmp, optmp); // (x*sqrt(3) - 1) / (x + sqrt(3))
        delete_UFs(optmp, tmp, NULL);
        res = calcAtan_UF(new_arg, new_arg, 1);
        Six = convertInteger_UF(6);
        optmp = div_UF(Pi, Six);
        res = call2_arg2(&add_UF,optmp,res);
        delete_UFs(optmp, new_arg, Six, NULL);
    }
    else 
        res = calcAtan_UF(arg, arg, 1);

    delete_UF(optmp2);

    optmp = div_UF(Pi, Two);
    if(compare_UF(res,optmp ) == 1)
    {
        delete_UF(res);
        res = div_UF(Pi, Two);
    }
    delete_UF(optmp);

    tmp = div_UF(Pi, Two);
    tmp->sign = (-1)*tmp->sign;
    if(compare_UF(tmp, res) == 1)
    { // if res<-pi/2 then return -pi/2
        delete_UF(res);
        res = div_UF(Pi, Two);
        res->sign = (-1)*res->sign;
    }

    delete_UFs(tmp, arg, eps, One, Two, SqrtThree, NULL);
    return res;
}

//////////////////////////////////////////////////////////////////////////
//                          atan2_UF                                    //
//////////////////////////////////////////////////////////////////////////
Unifloat* atan2_UF(Unifloat* x, Unifloat* y)
{
    /*
    res - contains current argument and may be returned as result
    pi,Zero, Two, Four - are constant numbers
    arg - is the same as x, using for x don't change in calculating
    eps - contains EPSILON
    tmp - auxiliary variable
    x_sign - sign of x
    */
    Unifloat* res, *tmp, *optmp, *optmp2;
    Unifloat* Zero;
    Unifloat* Two;
    Unifloat* Four;
    Unifloat* x_sign=convertInteger_UF(x->sign);

    res = clone(x);

    if(isZero_UF(x)&&(y->sign ==  - 1))
    {
        optmp = mul_UF(x_sign, Pi);
        delete_UFs(res, x_sign, NULL);
        return optmp;
    }

    if(isZero_UF(x)&&(y->sign == 1))
    {
        delete_UF(x_sign);
        return res;
    }

    if((x->sign ==  -1)&&isZero_UF(y))
    {
        optmp = convertInteger_UF(-2);
        optmp = call2_arg2(&div_UF, Pi, optmp);
        delete_UFs(res, x_sign, NULL);
        return optmp;
    }

    Two = convertInteger_UF(2);
    if((x->sign == 1)&&isZero_UF(y))
    {    
        optmp = div_UF(Pi, Two);
        delete_UFs(res, x_sign, Two, NULL);
        return optmp;
    }

    if(isNan_UF(x)||isNan_UF(y))
    {
        res->kind = NaN;
        delete_UFs(x_sign, Two, NULL);
        return res;
    }

    delete_UF(res);
   
    if (isZero_UF(x)&&isZero_UF(y)&&(y->sign==-1))
    {
        optmp = mul_UF(x_sign, Pi);
        delete_UFs(x_sign, Two, NULL);
        return optmp;
    }

    Zero = createZero_UF();
    if (isZero_UF(x)&&isZero_UF(y)&&(y->sign==1))
    {
        optmp = mul_UF(x_sign, Zero);
        delete_UFs(x_sign, Two, Zero, NULL);
        return optmp;
    }

    if((isInfinity_UF(x) != 1)
        &&(isInfinity_UF(y))&&(y->sign == -1))
    {
        optmp = mul_UF(x_sign, Pi);
        delete_UFs(x_sign, Two, Zero, NULL);
        return optmp;
    }

    if((isInfinity_UF(x) != 1)&&(isInfinity_UF(y))&&(y->sign == 1))
    {
        optmp =  mul_UF(x_sign, Zero);
        delete_UFs(x_sign, Two, Zero, NULL);
        return optmp;
    }

    delete_UF(Zero);
    if(isInfinity_UF(x)&&(isInfinity_UF(y) != 1))
    {
        optmp = mul_UF(Two, x_sign);
        
        optmp = call2_arg2(&div_UF, Pi, optmp);
        delete_UFs(x_sign, Two, NULL);
        return optmp;
    }

    Four = convertInteger_UF(4);
    if(isInfinity_UF(x)&&isInfinity_UF(y)&&(y->sign ==  -1))
    {
        tmp = convertInteger_UF(3*x->sign);
        optmp = div_UF(Pi, Four);
        optmp = call2_arg2(&mul_UF,tmp, optmp);
        delete_UFs(tmp, x_sign, Four, Two, NULL);
        return optmp;
    }
    
    if(isInfinity_UF(x)&&isInfinity_UF(y)&&(y->sign == 1))
    {
        optmp = mul_UF(Four,x_sign);
        optmp = call2_arg2(&div_UF, Pi, optmp);
        delete_UFs(x_sign, Four, Two, NULL);
        return optmp;
    }

    optmp = abs_UF(x);
    optmp2 = abs_UF(y);
    optmp = call2_arg1(&div_UF, optmp,  optmp2);
    res = atan_UF(optmp);

    if (y->sign == -1) res = call2_arg2(&sub_UF, Pi, res);
    if (x->sign == -1) res->sign = -res->sign;

    delete_UFs(optmp, optmp2, Four, x_sign, Two, NULL);
    return res;
}

//////////////////////////////////////////////////////////////////////////
//                          arrangeArgument_UF                          //
//////////////////////////////////////////////////////////////////////////
Unifloat* arrangeArgument_UF(Unifloat* x)
{
    Unifloat* int_part, *optmp, *optmp2;
    Unifloat* res;
    uint  i;
    Unifloat* Two = convertInteger_UF(2);
    Unifloat* TwoPi = mul_UF(Two, Pi);    
    res = clone(x);
    optmp2 = abs_UF(x);
    if(compare_UF(optmp2, TwoPi) == 1)
    {   
        int_part = div_UF(x, TwoPi);

        for(i = (int_part->exp) + 1; i <= PRECISION; i++)
            setMant_UF(int_part, i, 0);

        optmp = mul_UF(int_part, TwoPi);
        res = call2_arg1(&sub_UF,res, optmp);
        delete_UFs(optmp, int_part, NULL);
    }
    delete_UF(optmp2);    

    if(res->sign == -1) res = call2_arg1(&add_UF, res, TwoPi);

    delete_UFs(Two, TwoPi, NULL);
    return res;
}

//////////////////////////////////////////////////////////////////////////
//                          calcAtan_UF                                 //
//////////////////////////////////////////////////////////////////////////
Unifloat* calcAtan_UF(Unifloat* x, Unifloat* y, int i)
{
    Unifloat* tmp, * aux, * ind1, * ind2, *optmp;

    optmp = abs_UF(y);
    if(compareWithError_UF(optmp) != 1) 
    {
        delete_UF(optmp);
        return createZero_UF();
    }
    delete_UF(optmp);

    ind1 = convertInteger_UF(2 * i - 1);
    ind2 = convertInteger_UF(2 * i + 1);
    aux = mul_UF(x, x);

    optmp = mul_UF(y, aux);
    optmp = call2_arg1(&mul_UF, optmp, ind1);
    tmp = div_UF(optmp, ind2);
    delete_UF(optmp);

    tmp->sign = -tmp->sign;

    optmp = calcAtan_UF(x, tmp, i + 1);
    optmp = call2_arg2(&add_UF,y, optmp);
    
    delete_UFs(tmp, aux, ind1, ind2, NULL);
    return optmp;
}

//////////////////////////////////////////////////////////////////////////
//                          calcTan_UF                                  //
//////////////////////////////////////////////////////////////////////////
Unifloat* calcTan_UF(Unifloat* x, int i)
{
    Unifloat* tmp, *optmp, *optmp2, *res;

    if(i == ((int)(PRECISION/4) + 7)) 
        return createOne_UF();
    
    if(i == 1) 
    {
        optmp = calcTan_UF(x, i + 1);
        res = div_UF(x, optmp);
        delete_UF(optmp);
        return res;
    }
    else    
    {
        optmp = mul_UF(x, x);
        optmp2 = calcTan_UF(x, i + 1);
        tmp = div_UF(optmp, optmp2);
        delete_UFs(optmp, optmp2, NULL);

        optmp = convertInteger_UF(2*i - 3);
        res = sub_UF(optmp, tmp);
        delete_UFs(optmp, tmp, NULL);
        return res;
    }
}
