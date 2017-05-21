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

#include <string.h>
#include <limits.h>

#include <stdlib.h> /* NULL */
#include <stdarg.h> /* va_list */

#include "unifloat/unifloat.h"
#include "unifloat/unifloat_complex.h"
#include "unifloat/exp.h"
#include "unifloat/trig.h"
#include "unifloat/debug.h"

//////////////////////////////////////////////////////////////////////////
//                     create_UFComplex                                 //
//////////////////////////////////////////////////////////////////////////
UnifloatComplex * create_UFComplex( Unifloat* Re, Unifloat* Im)
{
    UnifloatComplex* new_UFComplex = (UnifloatComplex*)malloc(sizeof(UnifloatComplex));
    new_UFComplex->Re = Re;
    new_UFComplex->Im = Im;
    return new_UFComplex;
}

//////////////////////////////////////////////////////////////////////////
//                     delete_UFComplex                                 //
//////////////////////////////////////////////////////////////////////////
void delete_UFComplex( UnifloatComplex* x)
{
    delete_UF(x->Re);
    delete_UF(x->Im);
    free(x);
}

//////////////////////////////////////////////////////////////////////////
//                   delete_UFsComplex                                  //
//////////////////////////////////////////////////////////////////////////
void delete_UFsComplex(UnifloatComplex* p1, ...)
{
    UnifloatComplex* p;
    va_list ap;
    delete_UFComplex(p1);
    va_start(ap, p1);
    while((p = va_arg(ap, UnifloatComplex*)) != NULL) 
        delete_UFComplex(p);
    
    va_end(ap);
}

//////////////////////////////////////////////////////////////////////////
//                       clone_Complex                                  //
//////////////////////////////////////////////////////////////////////////
UnifloatComplex* clone_Complex(UnifloatComplex* x)
{
    UnifloatComplex* new_UFComplex = create_UFComplex(x->Re, x->Im);
    memcpy(new_UFComplex->Re, x->Re, sizeof(Unifloat));
    memcpy(new_UFComplex->Im, x->Im, sizeof(Unifloat));
    return new_UFComplex;
}

//////////////////////////////////////////////////////////////////////////
//                        copy_Complex                                  //
//////////////////////////////////////////////////////////////////////////
void copy_Complex(UnifloatComplex* src, UnifloatComplex* dst)
{
    memcpy(dst->Re, src->Re, sizeof(Unifloat));
    memcpy(dst->Im, src->Im, sizeof(Unifloat));
}

//////////////////////////////////////////////////////////////////////////
//                           abs_UFComplex                              //
//////////////////////////////////////////////////////////////////////////
Unifloat* abs_UFComplex(UnifloatComplex* x)
{
    Unifloat* res, *optmp;

    if (isNan_UF(x->Re) || isNan_UF(x->Im))
        return clone(nan_UF);

    if (isInfinity_UF(x->Re) || isInfinity_UF(x->Im))
        return clone(infinity_UF);

    res = mul_UF(x->Re, x->Re);

    optmp = mul_UF(x->Im, x->Im);
    res = call2_arg1(&add_UF, res, optmp);
    delete_UF(optmp);

    res = call1_arg1(&sqrt_UF, res);
    return res;
}

//////////////////////////////////////////////////////////////////////////
//                           carg_UF                                    //
//////////////////////////////////////////////////////////////////////////
Unifloat* carg_UF(UnifloatComplex* x) {
    return atan2_UF(x->Im, x->Re);
}

//////////////////////////////////////////////////////////////////////////
//                           add_UFComplex                              //
//////////////////////////////////////////////////////////////////////////
UnifloatComplex* add_UFComplex(UnifloatComplex* x, UnifloatComplex* y)
{
    UnifloatComplex* res;
    Unifloat* newRe;
    Unifloat* newIm;
    newRe = add_UF(x->Re, y->Re);
    newIm = add_UF(x->Im, y->Im);
    res = create_UFComplex(newRe, newIm);
    return res;
}

//////////////////////////////////////////////////////////////////////////
//                           sub_UFComplex                              //
//////////////////////////////////////////////////////////////////////////
UnifloatComplex* sub_UFComplex(UnifloatComplex* x, UnifloatComplex* y)
{
    UnifloatComplex* res;
    Unifloat* newRe;
    Unifloat* newIm;
    newRe = sub_UF(x->Re, y->Re);
    newIm = sub_UF(x->Im, y->Im);
    res = create_UFComplex(newRe, newIm);
    return res;
}

//////////////////////////////////////////////////////////////////////////
//                           mul_UFComplex                              //
//////////////////////////////////////////////////////////////////////////
UnifloatComplex* mul_UFComplex(UnifloatComplex* x, UnifloatComplex* y)
{
    UnifloatComplex* res;
    Unifloat* newRe;
    Unifloat* newIm;
    Unifloat* optmp, *optmp2;

    optmp = mul_UF(x->Re, y->Re);
    optmp2 = mul_UF(x->Im, y->Im);
    newRe = sub_UF(optmp, optmp2);
    delete_UFs(optmp, optmp2, NULL);

    optmp = mul_UF(x->Im, y->Re);
    optmp2 = mul_UF(x->Re, y->Im);
    newIm = add_UF(optmp, optmp2);
    delete_UFs(optmp, optmp2, NULL);

    res = create_UFComplex(newRe, newIm);
    return res;
}

//////////////////////////////////////////////////////////////////////////
//                           div_UFComplex                              //
//////////////////////////////////////////////////////////////////////////
UnifloatComplex* div_UFComplex(UnifloatComplex* x, UnifloatComplex* y)
{
    UnifloatComplex* res;
    Unifloat* newRe;
    Unifloat* newIm;
    Unifloat* optmp, *optmp2;

    optmp = mul_UF(x->Re, y->Re);
    optmp2 =  mul_UF(x->Im, y->Im);
    newRe = add_UF(optmp, optmp2);
    delete_UFs(optmp, optmp2, NULL);

    optmp = mul_UF(y->Re, y->Re);
    optmp2 = mul_UF(y->Im, y->Im);
    optmp = call2_arg1(&add_UF, optmp, optmp2);
    newRe = call2_arg1(&div_UF, newRe, optmp);
    delete_UFs(optmp, optmp2, NULL);

    optmp = mul_UF(x->Im, y->Re);
    optmp2 = mul_UF(x->Re, y->Im);
    newIm = sub_UF(optmp, optmp2);
    delete_UFs(optmp, optmp2, NULL);

    optmp = mul_UF(y->Re, y->Re);
    optmp2 = mul_UF(y->Im, y->Im);
    optmp = call2_arg1(&add_UF, optmp, optmp2);
    newIm = call2_arg1(&div_UF, newIm, optmp);
    res = create_UFComplex(newRe, newIm);
    return res;
}
