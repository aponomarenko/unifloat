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

#include "unifloat/unifloat.h"
#include "unifloat/exp.h"
#include "unifloat/cexp.h"
#include "unifloat/trig.h"
#include "unifloat/debug.h"

UnifloatComplex* cexp_UF(UnifloatComplex* x)
{
    UnifloatComplex* tmp1;
    Unifloat* exponent;
    UnifloatComplex* arg = clone_Complex(x);
    UnifloatComplex* res;
    Unifloat* Zero = createZero_UF();
    int  sign;

    if (x->Im->sign == -1)
    {
        x->Im->sign = 1;
        res = cexp_UF(x);
        x->Im->sign = -1;
        res->Im->sign = (-1) * res->Im->sign;
        return res;
    }    

    if (isZero_UF(x->Re) && isZero_UF(x->Im))
        return create_UFComplex(createOne_UF(), createZero_UF());

    if (isNormal_UF(x->Re) && isInfinity_UF(x->Im))
        return create_UFComplex(clone(nan_UF), clone(nan_UF));

    if (isNormal_UF(x->Re) && isNan_UF(x->Im))
        return create_UFComplex(clone(nan_UF), clone(nan_UF));

    if (isInfinity_UF(x->Re) && (x->Re->sign == 1) && isZero_UF(x->Im))
        return clone_Complex(x);

    if (isInfinity_UF(x->Re) && (x->Re->sign == -1) && isNormal_UF(x->Im))
        return create_UFComplex(mul_UF(Zero, cos_UF(x->Im)), mul_UF(Zero, sin_UF(x->Im)));

    if (isInfinity_UF(x->Re) && (x->Re->sign == 1) && isNormal_UF(x->Im) && !isZero_UF(x->Im))
        return create_UFComplex(mul_UF(infinity_UF, cos_UF(x->Im)), mul_UF(infinity_UF, sin_UF(x->Im)));
        
    if (isInfinity_UF(x->Re) && (x->Re->sign == -1) && isInfinity_UF(x->Im))
        return create_UFComplex(createZero_UF(), createZero_UF());

    if (isInfinity_UF(x->Re) && (x->Re->sign == 1) && isInfinity_UF(x->Im))
        return create_UFComplex(clone(infinity_UF), clone(nan_UF));

    if (isInfinity_UF(x->Re) && (x->Re->sign == -1) && isNan_UF(x->Im))
        return create_UFComplex(createZero_UF(), createZero_UF());

    if (isInfinity_UF(x->Re) && (x->Re->sign == 1) && isNan_UF(x->Im))
        return create_UFComplex(clone(infinity_UF), clone(nan_UF));

    if (isNan_UF(x->Re) && isZero_UF(x->Im))
        return clone_Complex(x);

    if (isNan_UF(x->Re) && !isNan_UF(x->Im) && !isZero_UF(x->Im))
        return create_UFComplex(clone(nan_UF), clone(nan_UF));

    if (isNan_UF(x->Re) && isNan_UF(x->Im))
        return create_UFComplex(clone(nan_UF), clone(nan_UF));

    exponent = exp_UF(arg->Re);

    tmp1 = create_UFComplex(cos_UF(arg->Im), sin_UF(arg->Im));

    res = clone_Complex(tmp1);
    res->Re = mul_UF(exponent, tmp1->Re);
    res->Im = mul_UF(exponent, tmp1->Im);

    if (isUnderflow_UF(res->Re))
    {
        sign = res->Re->sign;
        res->Re = createZero_UF();
        res->Re->sign = sign;
    }
    if (isUnderflow_UF(res->Im))
    {
        sign = res->Im->sign;
        res->Im = createZero_UF();
        res->Im->sign = sign;
    }
    
    if (isOverflow_UF(res->Re))
    {
        sign = res->Re->sign;
        res->Re = clone(infinity_UF);
        res->Re->sign = sign;
    }
    if (isOverflow_UF(res->Im))
    {
        sign = res->Im->sign;
        res->Im = clone(infinity_UF);
        res->Im->sign = sign;
    }

    return res;
}

UnifloatComplex* clog_UF(UnifloatComplex* x)
{
    UnifloatComplex* res = create_UFComplex(createZero_UF(),
                           createZero_UF());
    UnifloatComplex* arg = clone_Complex(x);
    Unifloat* PiFour = div_UF(Pi, convertInteger_UF(4));
    Unifloat* ThreePiFour = mul_UF(PiFour, convertInteger_UF(3));
    Unifloat* PiTwo = div_UF(Pi, convertInteger_UF(2));


    if (x->Im->sign == -1)
    {
        x->Im->sign = 1;
        res = clog_UF(x);
        x->Im->sign = -1;
        res->Im->sign = (-1) * res->Im->sign;
        return res;
    }    

    if (isZero_UF(x->Re) && (x->Re->sign == -1) && isZero_UF(x->Im))
    {
        res = create_UFComplex(clone(infinity_UF), clone(Pi));
        res->Re->sign = -1;
        return res;
    }

    if (isZero_UF(x->Re) && (x->Re->sign == 1) && isZero_UF(x->Im))
    {
        res = create_UFComplex(clone(infinity_UF), createZero_UF());
        res->Re->sign = -1;
        return res;
    }

    if (isNormal_UF(x->Re) && isInfinity_UF(x->Im))
        return create_UFComplex(clone(infinity_UF), clone(PiTwo));

    if (isNormal_UF(x->Re) && isNan_UF(x->Im))
        return create_UFComplex(clone(nan_UF), clone(nan_UF));

    if (isInfinity_UF(x->Re) && (x->Re->sign == -1) && isNormal_UF(x->Im))
        return create_UFComplex(clone(infinity_UF), clone(Pi));

    if (isInfinity_UF(x->Re) && (x->Re->sign == 1) && isNormal_UF(x->Im))
        return create_UFComplex(clone(infinity_UF), createZero_UF());

    if (isInfinity_UF(x->Re) && (x->Re->sign == -1) && isInfinity_UF(x->Im))
        return create_UFComplex(clone(infinity_UF), clone(ThreePiFour));

    if (isInfinity_UF(x->Re) && (x->Re->sign == 1) && isInfinity_UF(x->Im))
        return create_UFComplex(clone(infinity_UF), clone(PiFour));

    if (isInfinity_UF(x->Re) && isNan_UF(x->Im))
        return create_UFComplex(clone(infinity_UF), clone(nan_UF));

    if (isNan_UF(x->Re) && isNormal_UF(x->Im))
        return create_UFComplex(clone(nan_UF), clone(nan_UF));

    if (isNan_UF(x->Re) && isInfinity_UF(x->Im))
        return create_UFComplex(clone(infinity_UF), clone(nan_UF));

    if (isNan_UF(x->Re) && isNan_UF(x->Im))
        return create_UFComplex(clone(nan_UF), clone(nan_UF));

    res->Re = log_UF(abs_UFComplex(arg));
    res->Im = carg_UF(arg);

    return res;
}

UnifloatComplex* clog10_UF(UnifloatComplex* x)
{
    UnifloatComplex* res = create_UFComplex(createZero_UF(),
                           createZero_UF());
    UnifloatComplex* arg = clone_Complex(x);
    Unifloat* LN10 = log_UF(convertInteger_UF(10));
    Unifloat* PiLn10 = div_UF(Pi, LN10);
    Unifloat* PiLn10Four = div_UF(PiLn10, convertInteger_UF(4));
    Unifloat* ThreePiLn10Four = mul_UF(PiLn10Four, convertInteger_UF(3));
    Unifloat* PiLn10Two = div_UF(PiLn10, convertInteger_UF(2));
    
    if (x->Im->sign == -1)
    {
        x->Im->sign = 1;
        res = clog10_UF(x);
        x->Im->sign = -1;
        res->Im->sign = (-1) * res->Im->sign;
        return res;
    }    

    if (isZero_UF(x->Re) && (x->Re->sign == -1) && isZero_UF(x->Im))
    {
        res = create_UFComplex(clone(infinity_UF), clone(Pi));
        res->Re->sign = -1;
        return res;
    }

    if (isZero_UF(x->Re) && (x->Re->sign == 1) && isZero_UF(x->Im))
    {
        res = create_UFComplex(clone(infinity_UF), createZero_UF());
        res->Re->sign = -1;
        return res;
    }

    if (isNormal_UF(x->Re) && isInfinity_UF(x->Im))
        return create_UFComplex(clone(infinity_UF), clone(PiLn10Two));

    if (isNormal_UF(x->Re) && isNan_UF(x->Im))
        return create_UFComplex(clone(nan_UF), clone(nan_UF));

    if (isInfinity_UF(x->Re) && (x->Re->sign == -1) && isNormal_UF(x->Im))
        return create_UFComplex(clone(infinity_UF), clone(PiLn10));

    if (isInfinity_UF(x->Re) && (x->Re->sign == 1) && isNormal_UF(x->Im))
        return create_UFComplex(clone(infinity_UF), createZero_UF());

    if (isInfinity_UF(x->Re) && (x->Re->sign == -1) && isInfinity_UF(x->Im))
        return create_UFComplex(clone(infinity_UF), clone(ThreePiLn10Four));

    if (isInfinity_UF(x->Re) && (x->Re->sign == 1) && isInfinity_UF(x->Im))
        return create_UFComplex(clone(infinity_UF), clone(PiLn10Four));

    if (isInfinity_UF(x->Re) && isNan_UF(x->Im))
        return create_UFComplex(clone(infinity_UF), clone(nan_UF));

    if (isNan_UF(x->Re) && isNormal_UF(x->Im))
        return create_UFComplex(clone(nan_UF), clone(nan_UF));

    if (isNan_UF(x->Re) && isInfinity_UF(x->Im))
        return create_UFComplex(clone(infinity_UF), clone(nan_UF));

    if (isNan_UF(x->Re) && isNan_UF(x->Im))
        return create_UFComplex(clone(nan_UF), clone(nan_UF));

    res->Re = log_UF(abs_UFComplex(arg));
    res->Re = div_UF(res->Re, LN10);

    res->Im = carg_UF(arg);
    res->Im = div_UF(res->Im, LN10);
    
    return res;
}

UnifloatComplex* cpow_UF(UnifloatComplex* x, UnifloatComplex* y)
{
    UnifloatComplex* res = clone_Complex(x);
    UnifloatComplex* tmp1 = clone_Complex(x), *tmp2 = clone_Complex(x);

    if (isNan_UF(x->Re) && isInfinity_UF(x->Im) && 
        isNormal_UF(y->Re) && isNormal_UF(y->Re))
        return create_UFComplex(clone(infinity_UF), clone(nan_UF));

    if (isInfinity_UF(x->Re) && isNan_UF(x->Im) &&
        isNormal_UF(y->Re) && isNormal_UF(y->Im))
        return create_UFComplex(clone(infinity_UF), clone(nan_UF));


    tmp1 = clog_UF(x);
    tmp2 = mul_UFComplex(y, tmp1);
    res = cexp_UF(tmp2);

    return res;
}

UnifloatComplex* csqrt_UF(UnifloatComplex* x)
{
    UnifloatComplex* res = create_UFComplex(createZero_UF(), createZero_UF());
    UnifloatComplex* Half = create_UFComplex(div_UF(createOne_UF(), convertInteger_UF(2)), createZero_UF());

    if (x->Im->sign == -1)
    {
        x->Im->sign = 1;
        res = csqrt_UF(x);
        x->Im->sign = -1;
        res->Im->sign = (-1) * res->Im->sign;
        return res;
    }    

    if (isZero_UF(x->Re) && isZero_UF(x->Im))
        return create_UFComplex(createZero_UF(), createZero_UF());

    if (isInfinity_UF(x->Im))
        return create_UFComplex(clone(infinity_UF), clone(infinity_UF));

    if (isNormal_UF(x->Re) && isNan_UF(x->Im))
        return create_UFComplex(clone(nan_UF), clone(nan_UF));

    if (isInfinity_UF(x->Re) && (x->Re->sign == -1) && isNormal_UF(x->Im))
        return create_UFComplex(createZero_UF(), clone(infinity_UF));

    if (isInfinity_UF(x->Re) && (x->Re->sign == 1) && isNormal_UF(x->Im))
        return create_UFComplex(clone(infinity_UF), createZero_UF());

    if (isInfinity_UF(x->Re) && (x->Re->sign == -1) && isNan_UF(x->Im))
        return create_UFComplex(clone(nan_UF), clone(infinity_UF));

    if (isInfinity_UF(x->Re) && (x->Re->sign == 1) && isNan_UF(x->Im))
        return clone_Complex(x);

    if (isNan_UF(x->Re) && isNormal_UF(x->Im))
        return create_UFComplex(clone(nan_UF), clone(nan_UF));

    if (isNan_UF(x->Re) && isNan_UF(x->Im))
        return create_UFComplex(clone(nan_UF), clone(nan_UF));

    res = cpow_UF(x, Half);

    return res;
}
