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

#include <stdio.h> //perror
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include <stdlib.h> /* NULL */
#include <stdarg.h> /* va_list */

#include "unifloat/unifloat.h"
#include "unifloat/constants.h"
#include "unifloat/exp.h"
#include "unifloat/trig.h"

Unifloat* max_UF;
Unifloat* min_UF;
Unifloat* infinity_UF;
Unifloat* nan_UF;

Unifloat* Pi;
Unifloat* Gamma;
Unifloat* E;
Unifloat* Ln_2;
Unifloat* Log2_e;
Unifloat* Log10_e;
Unifloat* Log2_10;
Unifloat* Epsilon;

static int initialized = 0;
void initialize_UF()
{
    int i;
    CString* cstr;
    if( initialized == 1) return;
    initialized = 1;

    cstr = create_CString(PI);
    Pi = convertBinaryString_UF(cstr);
    delete_CString(cstr);
    
    cstr = create_CString(GAMMA);
    Gamma = convertString_UF(cstr);
    delete_CString(cstr);
    
    cstr = create_CString(_E);
    E = convertString_UF(cstr);
    delete_CString(cstr);
    
    cstr = create_CString(LN_2);
    Ln_2 = convertBinaryString_UF(cstr);
    delete_CString(cstr);
    
    cstr = create_CString(LOG2_E);
    Log2_e = convertBinaryString_UF(cstr);
    delete_CString(cstr);
    
    cstr = create_CString(LOG10_E);
    Log10_e = convertBinaryString_UF(cstr);
    delete_CString(cstr);
    
    cstr = create_CString(LOG2_10);
    Log2_10 = convertBinaryString_UF(cstr);
    delete_CString(cstr);
    
    cstr = create_CString(EPSILON);
    Epsilon = convertString_UF(cstr);
    delete_CString(cstr);

    infinity_UF = create_UF(1, 1, Infinity);
    nan_UF = create_UF(1, 1, NaN);

    max_UF = createOne_UF();
    for(i=0; i< MAX_SIZE_UNIFLOAT; i++)
        max_UF->mant[i] = UINT_MAX;

    min_UF = createOne_UF();
}

void finalize_UF()
{
    //finalizing math
    if( initialized == 0) return;
    initialized = 0;
    delete_UFs(Pi, Gamma, E, Ln_2, Log2_e, Log10_e, Log2_10, Epsilon, NULL);

    //finalizing unifloat internals
    delete_UFs(max_UF, min_UF, infinity_UF, nan_UF, NULL);
}

Unifloat* call1_arg1(caller_UF func, Unifloat* x)
{
    Unifloat* optmp = (*func) (x);
    delete_UF(x);
    return optmp;
}

Unifloat* call2_arg1(caller_UF_xy func, Unifloat* x, Unifloat* y)
{
    Unifloat* optmp = (*func) (x, y);
    delete_UF(x);
    return optmp;
}

Unifloat* call2_arg2(caller_UF_xy func, Unifloat* x, Unifloat* y)
{
    Unifloat* optmp = (*func) (x, y);
    delete_UF(y);
    return optmp;
}

#ifdef DEBUG_ON
long count_UF = 0;
#endif

//////////////////////////////////////////////////////////////////////////
//                      create_UF                                       //
//////////////////////////////////////////////////////////////////////////
Unifloat* create_UF(int sign, int exponent, UnifloatKind kind)
{   
    int i;
    Unifloat* new_UF = (Unifloat*)malloc(sizeof(Unifloat));
    if(new_UF == NULL)
    {
        perror("Malloc error on creating a Unifloat");
        exit(-1);
    }
    new_UF->sign = sign;
    new_UF->exp = exponent;
    for(i = 0; i < MAX_SIZE_UNIFLOAT; i++)
        new_UF->mant[i] = 0;

    new_UF->kind = kind;

    #ifdef DEBUG_ON
    count_UF++;
    #endif

    return new_UF;
}

//////////////////////////////////////////////////////////////////////////
//                      delete_UF                                       //
//////////////////////////////////////////////////////////////////////////
void delete_UF(Unifloat* u)
{   
    #ifdef DEBUG_ON
    count_UF--;
    #endif

    free(u);
}

//////////////////////////////////////////////////////////////////////////
//                      delete_UFs                                      //
//////////////////////////////////////////////////////////////////////////
void delete_UFs(Unifloat* p1, ...)
{
    Unifloat* p;
    va_list ap;
    delete_UF(p1);
    va_start(ap, p1);
    while((p = va_arg(ap, Unifloat*)) != NULL) 
        delete_UF(p);

    va_end(ap);
}

//////////////////////////////////////////////////////////////////////////
//                            clone                                     //
//////////////////////////////////////////////////////////////////////////
Unifloat* clone(Unifloat* x)
{
    int i;
    Unifloat* new_UF = create_UF(x->sign, x->exp, x->kind);
    for(i = 0; i < MAX_SIZE_UNIFLOAT; i++)
        new_UF->mant[i] = x->mant[i];

    return new_UF;
}

//////////////////////////////////////////////////////////////////////////
//                            copy                                      //
//////////////////////////////////////////////////////////////////////////
void copy(Unifloat* src, Unifloat* dst)
{
    int i;
    dst->sign = src->sign;
    dst->exp = src->exp;
    dst->kind = src->kind;

    for(i = 0; i < MAX_SIZE_UNIFLOAT; i++)
        dst->mant[i] = src->mant[i];
}

//////////////////////////////////////////////////////////////////////////
//                    isOverflow_UF                                     //
//////////////////////////////////////////////////////////////////////////
Bool isOverflow_UF(Unifloat* x)
{
    Unifloat* optmp;
    int r;
    if(x->kind == NaN) 
        return false;
    if(x->kind == Infinity) 
        return false;

    optmp = abs_UF(x);
    r = compare_UF(optmp, max_UF);
    delete_UF(optmp);
    return r == 1;
}

//////////////////////////////////////////////////////////////////////////
//                    isUnderflow_UF                                    //
//////////////////////////////////////////////////////////////////////////
Bool isUnderflow_UF(Unifloat* x)
{
    Unifloat* optmp;
    int r;
    if(x->kind == NaN)
        return false;
    if(isZero_UF(x) == true)
        return false;

    optmp = abs_UF(x);
    r = compare_UF(optmp, min_UF);
    delete_UF(optmp);
    return  r == -1;
}

//////////////////////////////////////////////////////////////////////////
//                    createZero_UF                                     //
//////////////////////////////////////////////////////////////////////////
Unifloat* createZero_UF() {
    return create_UF(1, 1, Normal);
}

//////////////////////////////////////////////////////////////////////////
//                      createOne_UF                                    //
//////////////////////////////////////////////////////////////////////////
Unifloat* createOne_UF()
{
    Unifloat* res;

    res = create_UF(1, 1, Normal);
    res->mant[0] = 0x1 << (sizeof(int) * 8 - 1);

    return res;
}

//////////////////////////////////////////////////////////////////////////
//                      compare_UF                                      //
//////////////////////////////////////////////////////////////////////////
int compare_UF(Unifloat* x, Unifloat* y)
{
    /*return  1 if x > y
     *return -1 if x < y
     *return  0 if x = y
     */

    int i, res;
    Unifloat* absx, * absy;

    if(isNan_UF(x) && isNan_UF(y))
        return 0;

    if((isNan_UF(x) && !isNan_UF(y))
     ||(!isNan_UF(x) && isNan_UF(y)))
        return 2;

    if(isZero_UF(x) && isZero_UF(y))
        return 0;
    if(isZero_UF(x))
        return -y->sign;
    if(isZero_UF(y)) 
        return x->sign;

    if(x->sign > y->sign) 
        return 1;
    else if(x->sign < y->sign) 
        return -1;

    if(x->sign == -1)
    {
        absx = abs_UF(x);
        absy = abs_UF(y);
        res = -compare_UF(absx, absy);
        delete_UFs(absx, absy, NULL);
        return res;
    }
    
    if(isInfinity_UF(x))
    {
        if(isInfinity_UF(y))
            return 0;
        else 
            return 1;
    }
    
    if(isInfinity_UF(y))
    {
        if(isInfinity_UF(x))
            return 0;
        else 
            return -1;
    }

    if(x->exp > y->exp) 
        return 1;
    else if(x->exp < y->exp) 
        return -1;
    else if(x->exp == y->exp)
    {
        for(i = 0; i < MAX_SIZE_UNIFLOAT; i++)
        {
            if(x->mant[i] > y->mant[i])
                return 1;
            else if(x->mant[i] < y->mant[i])
                return -1;
        }
    }
    return 0;
}

//////////////////////////////////////////////////////////////////////////
//                      isNormal_UF                                     //
//////////////////////////////////////////////////////////////////////////
Bool isNormal_UF(Unifloat* x)
{
    if(x->kind == Normal) 
        return true;
    else 
        return false;
}

//////////////////////////////////////////////////////////////////////////
//                      isInfinity_UF                                   //
//////////////////////////////////////////////////////////////////////////
Bool isInfinity_UF(Unifloat* x)
{
    if(x->kind == Infinity)
        return true;
    else 
        return false;
}

//////////////////////////////////////////////////////////////////////////
//                      isNan_UF                                        //
//////////////////////////////////////////////////////////////////////////
Bool isNan_UF(Unifloat* x)
{
    if(x->kind == NaN) 
        return true;
    else 
        return false;
}

//////////////////////////////////////////////////////////////////////////
//                      isZero_UF                                       //
//////////////////////////////////////////////////////////////////////////
Bool isZero_UF(Unifloat* x)
{
    int i;
    
    if(x->kind != Normal)
        return false;

    for(i = 0; i < MAX_SIZE_UNIFLOAT; i++)
        if(x->mant[i] != 0)
            return false;

    return true;
}

//////////////////////////////////////////////////////////////////////////
//                              abs_UF                                  //
//////////////////////////////////////////////////////////////////////////
Unifloat* abs_UF(Unifloat* x)
{
    Unifloat* y = clone(x);
    y->sign = 1;
    return y;
}

//////////////////////////////////////////////////////////////////////////
//                           convertUnifloat_Float                      //
//////////////////////////////////////////////////////////////////////////
float convertUnifloat_Float(Unifloat* x)
{ // TODO
    return (float)convertUnifloat_Double(x);
}

//////////////////////////////////////////////////////////////////////////
//                           convertUnifloat_Double                     //
//////////////////////////////////////////////////////////////////////////
double convertUnifloat_Double(Unifloat* x)
{
    double doubleInitial=0;
    int m;
    uint um;
    uint * bytes = (uint *)&doubleInitial;
    m = x->exp;
    m -=2;
    um = m;
    um = (um << 21) >> 1;
    if(x->exp<2)
        um &= ~(1<<(sizeof(int)*8-2));
    else
        um |= (1<<(sizeof(int)*8-2));

    bytes[1] &= 0xC00FFFFF;
    bytes[1] |= um;

    bytes[1] += (x->mant[0] << 1) >> 12;
    bytes[0] = (x->mant[0] << 21) + (x->mant[1] >> 11);

    if(bytes[0] == ~0x0u)
        bytes[1] += 0x1;

    if(x->sign<0)
        bytes[1] += (1<<(sizeof(int)*8 - 1));
    
    return doubleInitial;
}

//////////////////////////////////////////////////////////////////////////
//                           convertUnifloat_LongDouble                 //
//////////////////////////////////////////////////////////////////////////
long double convertUnifloat_LongDouble(Unifloat* x)
{ // TODO
    return (long double)convertUnifloat_Double(x);
}

//////////////////////////////////////////////////////////////////////////
//                      convertBinaryString_UF                          //
//////////////////////////////////////////////////////////////////////////
Unifloat* convertBinaryString_UF(CString* num)
{
    Unifloat* res = createZero_UF();
    CString* cstr, *number;
    int sign = 1, point_position = -1, i;
    char c;
    number = clone_CString(num);
    if (charAt_CString(number, 0) == '-')
    {
        sign = -1;
        cstr = substring_CString(number, 1, length_CString(number));
        delete_CString(number);
        number = cstr;
    }
    if (length_CString(number) > MAX_SIZE_UNIFLOAT*32)
    {
        cstr = substring_CString(number, 0, MAX_SIZE_UNIFLOAT* 32);
        delete_CString(number);
        number = cstr;
    }
    for (i = 0; (charAt_CString(number, i) != '1')&&(i < length_CString(number)); i++) 
        if (charAt_CString(number, i) == '.')
            point_position = i;

    if (point_position >= 0) 
    {
        res->exp = point_position - i + 1;
        point_position = -1;
    }
    else
        res->exp = length_CString(number);

    if (i > 0)
    {
        cstr = substring_CString(number, i, length_CString(number));
        delete_CString(number);
        number = cstr;
    }

    for (i = 0; i < length_CString(number); i++)
    {
        c = charAt_CString(number, i);
        if (c == '.')
            point_position = i;
        else if (c == '1')
        {
            if (point_position == -1)
                setMant_UF(res, i+1, 1);
            else
                setMant_UF(res, i, 1);
        }
    }
     
    if (point_position >= 0) 
        res->exp = point_position;
    
    res->sign = sign;
    delete_CString(number);
    return res;
}

//////////////////////////////////////////////////////////////////////////
//                           convertString_UF                           //
//////////////////////////////////////////////////////////////////////////
Unifloat* convertString_UF(CString* number)
{
    Unifloat* optmp;
    CString* cstr;
    Unifloat* digit[10];
    Unifloat* base;
    Unifloat* res = createZero_UF();
    Unifloat* power = createOne_UF();
    int i, exp = 0, err;
    char c;
    int point_position = 0;
    int sign = 1;
    digit[0] = createZero_UF();
    digit[1] = createOne_UF();
    for (i = 2; i < 10; i++)
        digit[i] = add_UF(digit[i - 1], digit[1]);

    base = add_UF(digit[9], digit[1]);

    if (charAt_CString(number, 0) == '-')
    {
        sign = -1;
        cstr= substring_CString(number, 1, length_CString(number));
        delete_CString(number);
        number = cstr;
    }
    if (charAt_CString(number, 0) == '+')
    {
        sign = 1;
        cstr = substring_CString(number, 1, length_CString(number));
        delete_CString(number);
        number = cstr;
    }
    if (indexOfChar_CString(number, 'e') != -1)
    {
        cstr = substring_CString(number, indexOfChar_CString(number, 'e') + 1, length_CString(number));
        optmp = convertString_UF(cstr);
        delete_CString(cstr);
            exp = convertUnifloat_Integer(optmp, &err);
        delete_UF(optmp);
            cstr = substring_CString(number, 0, indexOfChar_CString(number, 'e'));
        delete_CString(number);
        number = cstr;
    }

    // if(length_CString(number) > PRECISION / 3 + 1)
    //     number = substring_CString(number, 0, PRECISION / 3);

    for(i = length_CString(number) - 1; i >= 0; i--)
    {
        c = charAt_CString(number, i);

        if((c != '.') && (c != ','))
        {
            optmp = mul_UF(digit[c - 48], power);
            res = call2_arg1(&add_UF, res, optmp);
            delete_UF(optmp);
            power = call2_arg1(&mul_UF, power, base);
        }
        else
            point_position = length_CString(number) - i - 1;
    }
    delete_UF(power);
    power = createOne_UF();

    if (point_position - exp > 0)
    {
        for (i = point_position - exp; i > 0; i--)
            power = call2_arg1(&mul_UF, power, base);

        res = call2_arg1(&div_UF, res, power);
    }
    else
    {
        for (i = exp - point_position; i > 0; i--)
            power = call2_arg1(&mul_UF, power, base);

        res = call2_arg1(&mul_UF, res, power);
    }

    res->sign = sign;

    normalize_UF(res);
    delete_UFs(power, base, NULL);
    for (i=0; i<10; i++)
        delete_UF(digit[i]);
    return res;
}

//////////////////////////////////////////////////////////////////////////
//                           convertInteger_UF                          //
//////////////////////////////////////////////////////////////////////////
Unifloat* convertInteger_UF(int number)
{
    Unifloat* res = createZero_UF();
    
    if (number < 0)
    {
        res->sign = -1;
        res->mant[0] = ~number + 1;
    }
    else
        res->mant[0] = number;

    res->exp = sizeof(int) * 8;
    round_UF(res, sizeof(int) * 8 + 1);

    res = normalize_UF(res);
    return res;
}

//////////////////////////////////////////////////////////////////////////
//                           convertUnifloat_Integer                    //
//////////////////////////////////////////////////////////////////////////
int convertUnifloat_Integer(Unifloat* x, int* error)
{
    long res = 0;
    long power = 1;
    int i;
    
    *error = 0; // it means that it is possible to convert x to integer number

    if (isNan_UF(x) || (isInfinity_UF(x))) 
    {
        *error = 1;// it means that conversion is impossible
        return 0;
    }

    for (i = x->exp;i > 0;i--)
    {
        res += getMant_UF(x,i) * power;
        power = power*2;
    }

    if (x->exp > 0)
        i = x->exp+1;
    else
        i = 1;

    for (; i <= PRECISION; i++)
        if (getMant_UF(x, i) == 1)
            *error = 2; // it means that conversion isn't precisely

    if (res > INT_MAX)
        *error = 3; // it means that result isn't representable in 'int' type

    res = res * x->sign;
    return (int)res;
}

//////////////////////////////////////////////////////////////////////////
//                           convertFloat_UF                            //
//////////////////////////////////////////////////////////////////////////
Unifloat* convertFloat_UF(float number)
{ // TODO
    return convertDouble_UF((double)number);
}

//////////////////////////////////////////////////////////////////////////
//                           convertDouble_UF                           //
//////////////////////////////////////////////////////////////////////////
Unifloat* convertDouble_UF(double number)
{
    Unifloat* res;

    uint i, j;
    ulong *bytes, mask;
    uint bytesIndex;
    uint digExp, digMant;
    int maxExp;
    int exponent = 0x0;

    res = createZero_UF();

    digExp = digExp_DoubleT;
    digMant = digMant_DoubleT;
    maxExp = maxExp_DoubleT;

    bytes = (ulong *)&number;
    mask = 0x1;
    bytesIndex = 0;

    for (i = 0; i < digMant - 1; i++)
    {
        setMant_UF(res, digMant - i, (bytes[bytesIndex] & mask) / mask);
        
        mask <<= 1;
        if(mask == 0x0)
        {
            bytesIndex++;
            mask = 0x1;
        }
    }

    j = 0x1;
    for (i = 0; i < digExp; i++)
    {
        exponent += ((bytes[bytesIndex] & mask) / mask * j);
        
        mask <<= 1;
        if(mask == 0x0)
        {
            bytesIndex++;
            mask = 0x1;
        }
        j <<= 1;
    }

    if (exponent == maxExp * 2 - 1)   /* exponent == 0x111...11 */
    {
        setMant_UF(res, 1, 1);

        if (getMant_UF(res, 2) == 1)
            res->kind = NaN;
        else
            res->kind = Infinity;
    }
    else
    {
        res->kind = Normal;
        res->exp = exponent - maxExp + 2;

        if (exponent == 0x0)
        {
            setMant_UF(res, 1, 0);
            res->exp++;
            res = normalize_UF(res);
            res = round_UF(res, PRECISION);
        }
        else
            setMant_UF(res, 1, 1);
    }
    
    res->sign = 1 - 2 * ((bytes[bytesIndex] & mask) / mask);
    
    return res;
}

//////////////////////////////////////////////////////////////////////////
//                           convertLongDouble_UF                       //
//////////////////////////////////////////////////////////////////////////
Unifloat* convertLongDouble_UF(long double number)
{ // TODO
    return convertDouble_UF((double)number);
}

//////////////////////////////////////////////////////////////////////////
//                              round_UF                                //
//////////////////////////////////////////////////////////////////////////
Unifloat* round_UF(Unifloat* x, uint precision)
{
    uint localIndex = ((precision - 1) % (sizeof(int) * 8)) + 1;
    uint globalIndex = (precision - localIndex) / (sizeof(int) * 8);
    
    int i;
    Bool needRound = false;
    uint transposition;

    if (isNormal_UF(x) == false)
        return x;

    if (isZero_UF(x))
    {
        x->exp = 1;
        return x;
    }
    
    /* is rounding needed*/
    /* fulfilling with zeros*/
    if ((x->mant[globalIndex] << localIndex) != 0)
        needRound = true;

    for (i = globalIndex + 1; i < MAX_SIZE_UNIFLOAT; i++)
        if (x->mant[i] != 0)
            needRound = true; 

    /* rounding */
    if (needRound == false)
    {
        /* fulfilling with zeros*/
        x->mant[globalIndex] >>= sizeof(int) * 8 - localIndex;
        x->mant[globalIndex] <<= sizeof(int) * 8 - localIndex;

        for (i = globalIndex + 1; i < MAX_SIZE_UNIFLOAT; i++)
            x->mant[i] = 0x0;

        return x;
    }
    
    if (getMant_UF(x, precision + 1) == 1) 
        transposition = 1;
    else
        transposition = 0;

    if (transposition == 1)
    {
        transposition <<= (sizeof(int) * 8 - localIndex);

        for(i = globalIndex; i >= 0; i--)
        {
            if(x->mant[i] > UINT_MAX - transposition)
            {
                x->mant[i] += transposition;
                transposition = 1;
            }
            else
            {
                x->mant[i] += transposition;
                transposition = 0;
            }
        }
    }
    
    if (transposition == 1)
    {
        for (i = globalIndex; i >=0; i--)
        {
            x->mant[i] >>= 1;

            if (i > 0)
                x->mant[i] += x->mant[i - 1] << (sizeof(int) * 8 - 1);
            else
                x->mant[i] += 0x1 << (sizeof(int) * 8 - 1);
        }

        x->exp++;
    }
    
    /* fulfilling with zeros*/
    x->mant[globalIndex] >>= sizeof(int) * 8 - localIndex;
    x->mant[globalIndex] <<= sizeof(int) * 8 - localIndex;

    for (i = globalIndex + 1; i < MAX_SIZE_UNIFLOAT; i++)
        x->mant[i] = 0x0;

    return x;
}

//////////////////////////////////////////////////////////////////////////
//                          compareWithError_UF                         //
//////////////////////////////////////////////////////////////////////////
int compareWithError_UF(Unifloat* x) // compare with 10^(-1)*PRECISION
{
    int res;
    Unifloat* error = createOne_UF();

    error->exp = -150; //( - 1)*PRECISION;
    res = compare_UF(x, error);
    delete_UF(error);
    return res;
}

//////////////////////////////////////////////////////////////////////////
//                      compareWithPrecision_UF                         //
//////////////////////////////////////////////////////////////////////////
int compareWithPrecision_UF(Unifloat* x, Unifloat* y, int amount)
{
    int s;
    Unifloat* difference;
    amount -= COMPARE_PRECISION;

    if(isNan_UF(x))
    {
        if(isNan_UF(y))
            return 0;
        else
            return 2;
    }
    else if(isNan_UF(y))
        return 2;

    if(isInfinity_UF(x) && x->sign == 1)
    {
        if(isInfinity_UF(y) && y->sign == 1)
            return 0;
        else 
            return 1;
    }
    else if(isInfinity_UF(x) && x->sign == -1)
    {
        if(isInfinity_UF(y) && y->sign == -1)
            return 0;
        else 
            return -1;
    }
    else
    {
        if(isInfinity_UF(y) && y->sign == 1)
            return -1;
        else if(isInfinity_UF(y) && y->sign == -1)
            return 1;
    }

    if((x->exp > y->exp + 1) && !isZero_UF(x) && !isZero_UF(y))
        return x->sign;
    else if((x->exp < y->exp - 1) && !isZero_UF(x) && !isZero_UF(y))
        return -y->sign;
    else
    {
        if (isZero_UF(x))
        {
            if (y->exp <= -amount)
                return 0;
        }
        if (isZero_UF(y))
        {
            if (x->exp <= -amount)
                return 0;
        }
        difference = sub_UF(x, y);

        if(isZero_UF(difference)) 
        {
            delete_UF(difference);
            return 0;
        }
        if(difference->exp - x->exp - 1 <= - amount)
        {
            delete_UF(difference);
            return 0;
        }
        s = difference->sign;
        delete_UF(difference);
        return s;
    }
}

//////////////////////////////////////////////////////////////////////////
//                         setMant_UF                                   //
//////////////////////////////////////////////////////////////////////////
void setMant_UF(Unifloat* x, uint ind, uint bit)
{
    uint localIndex = ((ind - 1) % (sizeof(int) * 8)) + 1;
    uint globalIndex = (ind - localIndex) / (sizeof(int) * 8);
    uint temp = 0x1;

    temp <<= (sizeof(int) * 8 - localIndex);
    temp = ~temp;
    x->mant[globalIndex] &= temp;
    
    temp = bit;
    temp <<= (sizeof(int) * 8 - localIndex);
    x->mant[globalIndex] |= temp;
}

//////////////////////////////////////////////////////////////////////////
//                         getMant_UF                                   //
//////////////////////////////////////////////////////////////////////////

uint getMant_UF(Unifloat* x, uint ind)
{
    uint localIndex = ((ind - 1) % (sizeof(int) * 8)) + 1;
    uint globalIndex = (ind - localIndex) / (sizeof(int) * 8);
    uint res;

    res = x->mant[globalIndex];
    res >>= (sizeof(int) * 8 - localIndex);
    res &= 0x1;
    return res;
}

//////////////////////////////////////////////////////////////////////////
//                            changeSign_UF                             //
//////////////////////////////////////////////////////////////////////////
Unifloat* changeSign_UF(Unifloat* x, int s) 
{   
    /* set sign times s. If s = 0 then invert sign */
    Unifloat* res = clone(x);

    if(s == 0) 
        res->sign = -res->sign;
    else
        res->sign = res->sign * s;
    return res;
}

//////////////////////////////////////////////////////////////////////////
//                           normalize_UF                               //
//////////////////////////////////////////////////////////////////////////
Unifloat* normalize_UF(Unifloat* x)
{
    int i, j;
    int significantSize = (PRECISION - 1 - ((PRECISION - 1) % (sizeof(int) * 8))) / (sizeof(int) * 8) + 1;

    if(getMant_UF(x, 1) == 0 && !isZero_UF(x))
    {
        for(j = 0; x->mant[j] == 0x0; j++);

        if(j)
        {
            x->exp -= sizeof(int) * 8 * j;
            for(i = 0; i < significantSize + 1; i++)
                x->mant[i] = (i + j < MAX_SIZE_UNIFLOAT) ? x->mant[i + j] : 0x0;
        }

        for(j = 0; (x->mant[0] >> (sizeof(int) * 8 - 1 - j)) == 0; j++);

        if(j)
        {
            x->exp -= j;
            for(i = 0; i < significantSize + 1; i++)
            {
                x->mant[i] <<= j;
                x->mant[i] += x->mant[i + 1] >> (sizeof(int) * 8 - j);
            }
        }
    }
    
    return x;
}

//////////////////////////////////////////////////////////////////////////
//                              add_UF                                  //
//////////////////////////////////////////////////////////////////////////
Unifloat* add_UF(Unifloat* x, Unifloat* y)    // x + y
{
    int i;
    int significantSize = 
        (PRECISION - 1 - ((PRECISION - 1) % (sizeof(int) * 8))) / (sizeof(int) * 8) + 1;
    uint transposition;

    uint temp;
    Unifloat* res;

    int localDisplacement;
    int globalDisplacement;
    
    /* signs */
    if((x->sign == -1) && (y->sign == 1))
    {
        x->sign = 1;
        res = sub_UF(y, x);
        x->sign = -1;
        return res;
    }
    if((x->sign == 1) && (y->sign == -1))
    {
        y->sign = 1;
        res = sub_UF(x, y);
        y->sign = -1;
        return res;
    }

    /* exceptions */
    if(isNan_UF(x) || isNan_UF(y))              //if x or y is NaN
        return clone(nan_UF);

    if(isInfinity_UF(x) || isInfinity_UF(y))    //inf + inf = inf
    {
        res = clone(infinity_UF);
        res->sign = x->sign;
        return res;
    }

    if(isZero_UF(x))
        return clone(y);

    if(isZero_UF(y))
        return clone(x);
    
    transposition = 0;

    if(x->exp > y->exp)
    {
        res = clone(x);

        localDisplacement = (x->exp - y->exp) % (sizeof(int) * 8);
        globalDisplacement = ((x->exp - y->exp) - localDisplacement) / (sizeof(int) * 8);

        /* addition */
        for(i = significantSize; i >=0; i--)
        {
            temp = 0x0;

            if(i - globalDisplacement - 1 >= 0 && localDisplacement != 0)
                temp = y->mant[i - globalDisplacement - 1] << (sizeof(int) * 8 - localDisplacement);

            if(i - globalDisplacement >= 0)
                temp += y->mant[i - globalDisplacement] >> localDisplacement;

            res->mant[i] += temp + transposition;

            if(x->mant[i] > UINT_MAX - temp - transposition || (temp == UINT_MAX && transposition == 1))
                transposition = 1;
            else
                transposition = 0;
        }
    }
    else
    {
        res = clone(y);

        localDisplacement = (y->exp - x->exp) % (sizeof(int) * 8);
        globalDisplacement = ((y->exp - x->exp) - localDisplacement) / (sizeof(int) * 8);

        /* addition */
        for(i = significantSize; i >=0; i--)
        {
            temp = 0x0;

            if(i - globalDisplacement - 1 >= 0 && localDisplacement != 0)
                temp = x->mant[i - globalDisplacement - 1] << (sizeof(int) * 8 - localDisplacement);

            if(i - globalDisplacement >= 0)
                temp += x->mant[i - globalDisplacement] >> localDisplacement;

            res->mant[i] += temp + transposition;

            if(y->mant[i] > UINT_MAX - temp - transposition || (temp == UINT_MAX && transposition == 1))
                transposition = 1;
            else
                transposition = 0;
        }
    }

    /* normalization */
    if(transposition == 1)
    {
        for(i = significantSize; i >=0; i--)
        {
            res->mant[i + 1] += res->mant[i] << (sizeof(int) * 8 - 1);
            res->mant[i] >>= 1;
        }

        res->mant[0] += 0x1 << (sizeof(int) * 8 - 1);
        res->exp++;
    }

    res = round_UF(res, PRECISION);
    return res;
}

//////////////////////////////////////////////////////////////////////////
//                              sub_UF                                  //
//////////////////////////////////////////////////////////////////////////
Unifloat* sub_UF(Unifloat* x, Unifloat* y)
{
    int i, j;
    int significantSize = (PRECISION - 1 - ((PRECISION - 1) % (sizeof(int) * 8))) / (sizeof(int) * 8) + 1;
    
    uint transposition;

    uint temp;
    Unifloat* res;
    
    int localDisplacement;
    int globalDisplacement;

    /* signs */
    if((x->sign ==  -1) && (y->sign == 1))
    {
        y->sign = -1;
        res = add_UF(x, y);
        y->sign = 1;
        return res;
    }
    if((x->sign == 1) && (y->sign == -1))
    {
        y->sign = 1;
        res = add_UF(x, y);
        y->sign = -1;
        return res;
    }

    /* exceptions */

    /* if x or y - NaN */
    if(isNan_UF(x) || isNan_UF(y))
        return clone(nan_UF);

    /* +- (inf - inf) = NaN */
    if(isInfinity_UF(x) &&  isInfinity_UF(y))
        return clone(nan_UF);

    /* +- (inf - norm) =  +- (inf) */
    if(isInfinity_UF(x) &&  isNormal_UF(y))
    {
        res = clone(infinity_UF);
        res->sign = x->sign;
        return res;
    }

    /* +- (norm - inf) = -+ (inf) */
    if(isInfinity_UF(y) &&  isNormal_UF(x))
    {
        res = clone(infinity_UF);
        res->sign = -y->sign;
        return res;
    }

    /* +- (|x| - |y|);|x| >= |y| */
    if(compare_UF(x, y) * x->sign > -1) 
    {
        res = clone(x);

        if(isZero_UF(y))
            return res;

        localDisplacement = (x->exp - y->exp) % (sizeof(int) * 8);
        globalDisplacement = ((x->exp - y->exp) - localDisplacement) / (sizeof(int) * 8);

        transposition = 0;

        /* subtraction */
        for(i = significantSize; i >= 0; i--)
        {
            temp = 0x0;

            if(i - globalDisplacement - 1 >= 0 && localDisplacement != 0)
                temp = y->mant[i - globalDisplacement - 1] << (sizeof(int) * 8 - localDisplacement);

            if(i - globalDisplacement >= 0)
                temp += y->mant[i - globalDisplacement] >> localDisplacement;

            res->mant[i] -= (temp + transposition);

            if((x->mant[i] < temp + transposition) || (temp == UINT_MAX && transposition == 1))
                transposition = 1;
            else
                transposition = 0;
        }
    }
    else
    {
        res = clone(y);
        res->sign = -res->sign;

        if(isZero_UF(x))
            return res;

        localDisplacement = (y->exp - x->exp) % (sizeof(int) * 8);
        globalDisplacement = ((y->exp - x->exp) - localDisplacement) / (sizeof(int) * 8);

        transposition = 0;

        /* subtraction */
        for(i = significantSize; i >= 0; i--)
        {
            temp = 0x0;

            if(i - globalDisplacement - 1 >= 0 && localDisplacement != 0)
                temp = x->mant[i - globalDisplacement - 1] << (sizeof(int) * 8 - localDisplacement);

            if(i - globalDisplacement >= 0)
                temp += x->mant[i - globalDisplacement] >> localDisplacement;

            res->mant[i] -= (temp + transposition);

            if((y->mant[i] < temp + transposition) || (temp == UINT_MAX && transposition == 1))
                transposition = 1;
            else
                transposition = 0;
        }
    }

    /* normalization */
    if(getMant_UF(res, 1) == 0 && !isZero_UF(res))
    {
        for(j = 0; res->mant[j] == 0x0; j++);

        if(j)
        {
            res->exp -= sizeof(int) * 8 * j;
            for(i = 0; i < significantSize + 1; i++)
                res->mant[i] = (i + j < MAX_SIZE_UNIFLOAT) ? res->mant[i + j] : 0x0;
        }

        for(j = 0; (res->mant[0] >> (sizeof(int) * 8 - 1 - j)) == 0; j++);

        if(j)
        {
            res->exp -= j;
            for(i = 0; i < significantSize + 1; i++)
            {
                res->mant[i] <<= j;
                res->mant[i] += res->mant[i + 1] >> (sizeof(int) * 8 - j);
            }
        }
    }

    res = round_UF(res, PRECISION);
    return res;
}

//////////////////////////////////////////////////////////////////////////
//                              mul_UF                                  //
//////////////////////////////////////////////////////////////////////////
Unifloat* mul_UF(Unifloat* x, Unifloat* y)
{
    int i, j, k;
    int significantSize = (PRECISION - 1 - ((PRECISION - 1) % (sizeof(int) * 8))) / (sizeof(int) * 8) + 1;

    uint temp;
    uint low = 0x0;
    Unifloat* res;
    /* if x or y - NaN */
    if(isNan_UF(x) || isNan_UF(y))
        return clone(nan_UF);

    /* if x and y - infinity then return infinity */
    if(isInfinity_UF(x))
    {
        if(isZero_UF(y))
            res = clone(nan_UF);
        else
        {
            res = clone(infinity_UF);
            res->sign = x->sign * y->sign;
            return res;
        }    
    }
    
    if(isInfinity_UF(y))
    {
        if(isZero_UF(x))
            res = clone(nan_UF);
        else
        {
            res = clone(infinity_UF);
            res->sign = x->sign * y->sign;
            return res;
        }
    }   
    res = createZero_UF();
    if(isZero_UF(x) || isZero_UF(y))
        return res;

    res->sign = x->sign * y->sign;
    res->exp = x->exp + y->exp;

    for(i = significantSize * 2 - 1; i >=0; i--)
    {
        for(j = significantSize * 2 - 1; j >=0; j--)
        {
            temp = 0x0;
            low = 0x0;
            low = ~low;
            low >>= (sizeof(int) * 4);

            temp = ((i % 2) ? (x->mant[(i - 1) / 2] & low) : (x->mant[i / 2] >> (sizeof(int) * 4))) 
                * ((j % 2) ? (y->mant[(j - 1) / 2] & low) : (y->mant[j / 2] >> (sizeof(int) * 4)));
            
            if((i + j) % 2 == 0)
            {
                k = (i + j) / 2;
                while(res->mant[k] > UINT_MAX - temp)
                {
                    res->mant[k] += temp;
                    temp = 0x1;
                    k--;
                }
                res->mant[k] += temp;
            }
            else
            {
                low = (temp & low) << (sizeof(int) * 4);

                k = (i + j + 1) / 2;
                while(res->mant[k] > UINT_MAX - low)
                {
                    res->mant[k] += low;
                    low = 0x1;
                    k--;
                }
                res->mant[k] += low;

                temp >>= (sizeof(int) * 4);

                k = (i + j - 1) / 2;
                while(res->mant[k] > UINT_MAX - temp)
                {
                    res->mant[k] += temp;
                    temp = 0x1;
                    k--;
                }
                res->mant[k] += temp;
            }
        }
    }
    
    /* normalization */
    if(getMant_UF(res, 1) == 0)
    {
        for(i = 0; i < significantSize + 1; i++)
        {
            res->mant[i] <<= 1;
            res->mant[i] += res->mant[i + 1] >> (sizeof(int) * 8 - 1);
        }

        res->exp--;
    }
    
    res = round_UF(res, PRECISION);
    return res;
}

//////////////////////////////////////////////////////////////////////////
//                                div_UF                                //
//////////////////////////////////////////////////////////////////////////
Unifloat* div_UF(Unifloat* x, Unifloat* y)
{
    Unifloat *res, *two, *initial, *denom, *temp, *optmp;

    double doubleInitial = 1.0;
    uint *bytes;
    int inaccuracy;

    if(isNan_UF(x) || isNan_UF(y))
        return clone(nan_UF);

    if(isInfinity_UF(x) && isInfinity_UF(y))
        return clone(nan_UF);

    if(isInfinity_UF(x) && isNormal_UF(y) )
     {
        temp = clone(infinity_UF);
        temp->sign = x->sign * y->sign;
        return temp;
     }

    if(isNormal_UF(x) && isInfinity_UF(y))
        return createZero_UF();

    if(isZero_UF(x) && isZero_UF(y))
        return clone(nan_UF);
    
    if(!isZero_UF(x) && isZero_UF(y))
     {
        temp = clone(infinity_UF);
        temp->sign = x->sign * y->sign;
        return temp;
     }

    if(isZero_UF(x))
        return createZero_UF();

    initial = createOne_UF();
    bytes = (uint *)&doubleInitial;

    bytes[1] += (y->mant[0] << 1) >> 12;
    bytes[0] = (y->mant[0] << 21) + (y->mant[1] >> 11);

    if(bytes[0] == ~0x0u)
        bytes[1] += 0x1;
    bytes[0] += 0x1;


    doubleInitial = 1.0 / doubleInitial;
    initial->mant[0] += ((bytes[1] << 12) >> 1) + (bytes[0] >> 21);
    initial->mant[1] = bytes[0] << 11;
    initial->exp = ((bytes[1] << 1) >> 21) - maxExp_DoubleT + 2;
    
    two = createOne_UF();
    two->exp++;

    denom = clone(y);
    denom->sign = 1;
    denom->exp = 1;

    for(inaccuracy = digMant_DoubleT; inaccuracy < PRECISION; inaccuracy *= 2)
    {
        optmp = mul_UF(initial, denom);
        temp = sub_UF(two, optmp);
        delete_UF(optmp);
        initial = call2_arg2(&mul_UF, temp, initial);
        delete_UF(temp);
    }

    res = mul_UF(initial, x);
    res->sign *= y->sign;
    res->exp -= y->exp - 1;

    delete_UFs(initial, two, denom, NULL);
    return res;
}

//////////////////////////////////////////////////////////////////////
//                          power_UF                                //
//////////////////////////////////////////////////////////////////////
Unifloat* power_UF(Unifloat* x, int n)
{
    int i;
    Unifloat* res, * tmp;
    res = convertInteger_UF(1);
    for(i = 1; i < n + 1; i++)
    {
        tmp = mul_UF(res, x);
        delete_UF(res);
        res = tmp;
    }

    return res;
}

//////////////////////////////////////////////////////////////////////
//                          factorial_UF                            //
//////////////////////////////////////////////////////////////////////
Unifloat* factorial_UF(int n)
{
    int i, res = 1;

    for(i = 1; i < n + 1; i++)
        res *= i;

    return convertInteger_UF(res);
}

//////////////////////////////////////////////////////////////////////
//                          print_UF                                //
//////////////////////////////////////////////////////////////////////
void print_UF(Unifloat* x) {
    printf("%.17g", convertUnifloat_Double(x));
}

//////////////////////////////////////////////////////////////////////////
//                                call_UF                                //
//////////////////////////////////////////////////////////////////////////
double call_UF(caller_UF func, double x)
{
    Unifloat* ux;
    Unifloat* ures;
    double res;

    ux = convertDouble_UF(x);
    ures = (*func)(ux);
    res = convertUnifloat_Double(ures);
    delete_UFs(ux, ures, NULL);

    return res;
}

//////////////////////////////////////////////////////////////////////////
//                                call_UF_nx                            //
//////////////////////////////////////////////////////////////////////////
double call_UF_nx(caller_UF_nx func, int n, double x)
{
    Unifloat* ux;
    Unifloat* ures;
    double res;

    ux = convertDouble_UF(x);
    ures = (*func)(n, ux);
    res = convertUnifloat_Double(ures);
    delete_UFs(ux, ures, NULL);

    return res;
}
