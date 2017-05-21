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
#include "unifloat/trig.h"
#include "unifloat/ctrig.h"
#include "unifloat/exp.h"
#include "unifloat/debug.h"

UnifloatComplex* cacosh_UF(UnifloatComplex* x)
{
    UnifloatComplex* aux1, * aux2;
    UnifloatComplex* res;
    Unifloat* real, * imag, * arg1, * arg2, * arg, * RE, * IM;
    Unifloat* One = createOne_UF();
    Unifloat* PiFour = div_UF(Pi, convertInteger_UF(4));
    Unifloat* ThreePiFour = mul_UF(PiFour, convertInteger_UF(3));
    Unifloat* PiTwo = div_UF(Pi, convertInteger_UF(2));

    if (x->Im->sign == -1) 
    {
        x->Im->sign = 1;
        res = cacosh_UF(x);
        x->Im->sign = -1;
        res->Im->sign = (-1) * res->Im->sign;
        return res;
    }
    
    if (isZero_UF(x->Re) && isZero_UF(x->Im))
        return create_UFComplex(createZero_UF(), clone(PiTwo));

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

    if ((compare_UF(x->Re, One) == 0) && isZero_UF(x->Im))
        return create_UFComplex(createZero_UF(), createZero_UF());

    if ((compare_UF(abs_UF(x->Re), One) == 0) && isZero_UF(x->Im))    
        return create_UFComplex(createZero_UF(), clone(Pi));
    
    aux1 = clone_Complex(x);
    aux1->Re = add_UF(aux1->Re, One);
    arg1 = acos_UF(div_UF(aux1->Re, abs_UFComplex(aux1)));

    aux2 = clone_Complex(x);
    aux2->Re = sub_UF(aux2->Re, One);
    arg2 = acos_UF(div_UF(aux2->Re, abs_UFComplex(aux2)));

    arg = add_UF(arg1, arg2);
    arg->exp--;

    real = sqrt_UF(abs_UFComplex(aux1));
    real = mul_UF(real, sqrt_UF(abs_UFComplex(aux2)));
    imag = clone(real);
    real = mul_UF(real, cos_UF(arg));
    imag = mul_UF(imag, sin_UF(arg));
    real = add_UF(x->Re, real);
    imag = add_UF(x->Im, imag);

    RE = mul_UF(real, real);
    RE = add_UF(RE, mul_UF(imag, imag));
    RE = sqrt_UF(RE);
    RE = log_UF(RE);

    IM = mul_UF(real, real);
    IM = add_UF(IM, mul_UF(imag, imag));
    IM = sqrt_UF(IM);
    IM = div_UF(real, IM);
    IM = acos_UF(IM);

//    aux1 = add_UFComplex(x, OneC);
//    aux2 = sub_UFComplex(OneC, x);
//
//    res = clog_UF(div_UFComplex(aux1, aux2));
//    res->Re = div_UF(res->Re, Two);
//    res->Im = div_UF(res->Im, Two);

    return create_UFComplex(RE, IM);
}

UnifloatComplex* casinh_UF(UnifloatComplex* x)
{
    UnifloatComplex* aux, * OneC;
    Unifloat* real, * imag, * temp, * RE, * IM;
    Unifloat* Zero = createZero_UF();
    Unifloat* One = createOne_UF();
    UnifloatComplex* res;
    Unifloat* PiFour = div_UF(Pi, convertInteger_UF(4));
    Unifloat* PiTwo = div_UF(Pi, convertInteger_UF(2));

    if (x->Re->sign == -1)
    {
        x->Re->sign = 1;
        x->Im->sign = (-1) * x->Im->sign;
        res = casinh_UF(x);
        x->Re->sign = -1;
        x->Im->sign = (-1) * x->Im->sign;
        res->Im->sign = (-1) * res->Im->sign;
        res->Re->sign = (-1) * res->Re->sign;
        return res;
    }

    if (x->Im->sign == -1)
    {
        x->Im->sign = 1;
        res = casinh_UF(x);
        x->Im->sign = -1;
        res->Im->sign = (-1) * res->Im->sign;
        return res;
    }
    if (isZero_UF(x->Re) && isZero_UF(x->Im))
        return create_UFComplex(createZero_UF(), createZero_UF());

    if (isNormal_UF(x->Re) && isInfinity_UF(x->Im))
        return create_UFComplex(clone(infinity_UF), clone(PiTwo));

    if (isNormal_UF(x->Re) && isNan_UF(x->Im))
        return create_UFComplex(clone(nan_UF), clone(nan_UF));

    if (isInfinity_UF(x->Re) && isNormal_UF(x->Im))
        return create_UFComplex(clone(infinity_UF), createZero_UF());

    if (isInfinity_UF(x->Re) && isInfinity_UF(x->Im))
        return create_UFComplex(clone(infinity_UF), clone(PiFour));

    if (isInfinity_UF(x->Re) && isNan_UF(x->Im))
        return create_UFComplex(clone(infinity_UF), clone(nan_UF));

    if (isNan_UF(x->Re) && isZero_UF(x->Im))
        return create_UFComplex(clone(nan_UF), createZero_UF());

    if (isNan_UF(x->Re) && isNormal_UF(x->Im) && !isZero_UF(x->Im))
        return create_UFComplex(clone(nan_UF), clone(nan_UF));

    if (isNan_UF(x->Re) && isInfinity_UF(x->Im))
        return create_UFComplex(clone(infinity_UF), clone(nan_UF));

    if (isNan_UF(x->Re) && isNan_UF(x->Im))
        return create_UFComplex(clone(nan_UF), clone(nan_UF));

    if (isZero_UF(x->Re) && (compare_UF(x->Im, One) == 0))
        return create_UFComplex(createZero_UF(), clone(PiTwo));

    OneC = create_UFComplex(One, Zero);
    aux = mul_UFComplex(x, x);
    aux = add_UFComplex(aux, OneC);
    
    real = abs_UFComplex(aux);
    temp = acos_UF(div_UF(aux->Re, real));
    temp->exp--;
    real = sqrt_UF(real);
    imag = clone(real);
    real = mul_UF(real, cos_UF(temp));
    imag = mul_UF(imag, sin_UF(temp));
    real = add_UF(x->Re, real);
    imag =add_UF(x->Im, imag);
    
    RE = mul_UF(real, real);
    RE = add_UF(RE, mul_UF(imag, imag));
    RE = sqrt_UF(RE);
    RE = log_UF(RE);

    IM = mul_UF(real, real);
    IM = add_UF(IM, mul_UF(imag, imag));
    IM = sqrt_UF(IM);
    IM = div_UF(real, IM);
    IM = acos_UF(IM);

    return create_UFComplex(RE, IM);
}

UnifloatComplex* catanh_UF(UnifloatComplex* x)
{
    UnifloatComplex *tempC, *OneC;
    Unifloat* real, * imag, * tmp;
    Unifloat* Zero = createZero_UF();
    Unifloat* One = createOne_UF();
    UnifloatComplex* res;
    Unifloat* PiTwo = div_UF(Pi, convertInteger_UF(2));

    if (x->Re->sign == -1)
    {
        x->Re->sign = 1;
        x->Im->sign = (-1) * x->Im->sign;
        res = catanh_UF(x);
        x->Re->sign = -1;
        x->Im->sign = (-1) * x->Im->sign;
        res->Im->sign = (-1) * res->Im->sign;
        res->Re->sign = (-1) * res->Re->sign;
        return res;
    }

    if (x->Im->sign == -1)
    {
        x->Im->sign = 1;
        res = catanh_UF(x);
        x->Im->sign = -1;
        res->Im->sign = (-1) * res->Im->sign;
        return res;
    }

    if (isZero_UF(x->Re) && isZero_UF(x->Im))
        return create_UFComplex(createZero_UF(), createZero_UF());

    if (isZero_UF(x->Re) && isNan_UF(x->Im))
        return create_UFComplex(createZero_UF(), clone(nan_UF));

    if ((compare_UF(x->Re, One) == 0) && isZero_UF(x->Im))
        return create_UFComplex(clone(infinity_UF), createZero_UF());

    if (isNormal_UF(x->Re) && isInfinity_UF(x->Im))
        return create_UFComplex(createZero_UF(), clone(PiTwo));

    if (isNormal_UF(x->Re) && !isZero_UF(x->Re) && isNan_UF(x->Im))
        return create_UFComplex(clone(nan_UF), clone(nan_UF));

    if (isInfinity_UF(x->Re) && isNormal_UF(x->Im))
        return create_UFComplex(createZero_UF(), clone(PiTwo));

    if (isInfinity_UF(x->Re) && isInfinity_UF(x->Im))
        return create_UFComplex(createZero_UF(), clone(PiTwo));

    if (isInfinity_UF(x->Re) && isNan_UF(x->Im))
        return create_UFComplex(createZero_UF(), clone(nan_UF));

    if (isNan_UF(x->Re) && isNormal_UF(x->Im))
        return create_UFComplex(clone(nan_UF), clone(nan_UF));

    if (isNan_UF(x->Re) && isInfinity_UF(x->Im))
        return create_UFComplex(createZero_UF(), clone(PiTwo));

    if (isNan_UF(x->Re) && isNan_UF(x->Im))
        return create_UFComplex(clone(nan_UF), clone(nan_UF));

    tempC = create_UFComplex(Zero, Zero);
    OneC = create_UFComplex(One, Zero);
    tempC = div_UFComplex(add_UFComplex(OneC, x), 
                                sub_UFComplex(OneC, x));
    
    real = log_UF(abs_UFComplex(tempC));
    
    tmp = sqrt_UF(add_UF(mul_UF(tempC->Im, tempC->Im), 
          mul_UF(tempC->Re, tempC->Re)));
    imag = acos_UF(div_UF(tempC->Re, tmp));
    
    real->exp--;
    imag->exp--;

    return create_UFComplex(real, imag);
}

UnifloatComplex* ccosh_UF(UnifloatComplex* x)
{
    Unifloat* real, * imag, * One;
    UnifloatComplex* res;
    One = createOne_UF();

    if (x->Re->sign == -1)
    {
        x->Re->sign = 1;
        x->Im->sign = (-1) * x->Im->sign;
        res = ccosh_UF(x);
        x->Re->sign = -1;
        x->Im->sign = (-1) * x->Im->sign;
        return res;
    }

    if (x->Im->sign == -1)
    {
        x->Im->sign = 1;
        res = ccosh_UF(x);
        x->Im->sign = -1;
        res->Im->sign = (-1)* res->Im->sign;
        return res;
    }

    if (isZero_UF(x->Re) && isZero_UF(x->Im))
        return create_UFComplex(createOne_UF(), createZero_UF());

    if (isZero_UF(x->Re) && isInfinity_UF(x->Im))
        return create_UFComplex(clone(nan_UF), createZero_UF());

    if (isZero_UF(x->Re) && isNan_UF(x->Im))
        return create_UFComplex(clone(nan_UF), createZero_UF());

    if (isNormal_UF(x->Re) && !isZero_UF(x->Re) && isInfinity_UF(x->Im))
        return create_UFComplex(clone(nan_UF), clone(nan_UF));

    if (isNormal_UF(x->Re) && !isZero_UF(x->Re) && isNan_UF(x->Im))
        return create_UFComplex(clone(nan_UF), clone(nan_UF));

    if (isInfinity_UF(x->Re) && isZero_UF(x->Im))
        return clone_Complex(x);

    if (isInfinity_UF(x->Re) && isNormal_UF(x->Im)&& !isZero_UF(x->Im))
        return create_UFComplex(mul_UF(infinity_UF, cos_UF(x->Im)),
                                     mul_UF(infinity_UF, sin_UF(x->Im)));

    if (isInfinity_UF(x->Re) && isInfinity_UF(x->Im))
        return create_UFComplex(clone(infinity_UF), clone(nan_UF));

    if (isInfinity_UF(x->Re) && isNan_UF(x->Im))
        return create_UFComplex(clone(infinity_UF), clone(nan_UF));

    if (isNan_UF(x->Re) && isZero_UF(x->Im))
        return clone_Complex(x);

    if (isNan_UF(x->Re) && !isZero_UF(x->Im))
        return create_UFComplex(clone(nan_UF), clone(nan_UF));

    if (isNan_UF(x->Re) && isNan_UF(x->Im))
        return create_UFComplex(clone(nan_UF), clone(nan_UF));

    real = exp_UF(x->Re);
    real = add_UF(real, div_UF(One, real));
    real = mul_UF(real, cos_UF(x->Im));
    real->exp--;

    imag = exp_UF(x->Re);
    imag = sub_UF(imag, div_UF(One, imag));
    imag = mul_UF(imag, sin_UF(x->Im));
    imag->exp--;
    res = create_UFComplex(real, imag);

    return res;
}

UnifloatComplex* csinh_UF(UnifloatComplex* x)
{
    Unifloat* real, * imag, * One;
    UnifloatComplex* res;
    
    One = createOne_UF();

    if (isInfinity_UF(x->Re) && isInfinity_UF(x->Im))
        return create_UFComplex(clone(infinity_UF), clone(nan_UF));

    if (isInfinity_UF(x->Re) && isNan_UF(x->Im))
        return create_UFComplex(clone(infinity_UF), clone(nan_UF));

    if (x->Re->sign == -1)
    {
        x->Re->sign = 1;
        x->Im->sign = (-1) * x->Im->sign;
        res = csinh_UF(x);
        x->Re->sign = -1;
        x->Im->sign = (-1) * x->Im->sign;
        res->Im->sign = (-1) * res->Im->sign;
        res->Re->sign = (-1) * res->Re->sign;
        return res;
    }

    if (x->Im->sign == -1)
    {
        x->Im->sign = 1;
        res = csinh_UF(x);
        x->Im->sign = -1;
        res->Im->sign = (-1) * res->Im->sign;
        return res;
    }

    if (isZero_UF(x->Re) && isZero_UF(x->Im))
        return clone_Complex(x);

    if (isZero_UF(x->Re) && isInfinity_UF(x->Im))
        return create_UFComplex(createZero_UF(), clone(nan_UF));

    if (isZero_UF(x->Re) && isNan_UF(x->Im))
        return create_UFComplex(createZero_UF(), clone(nan_UF));

    if (isNormal_UF(x->Re) && isInfinity_UF(x->Im))
        return create_UFComplex(clone(nan_UF), clone(nan_UF));

    if (isNormal_UF(x->Re) && !isZero_UF(x->Re) && isNan_UF(x->Im))
        return create_UFComplex(clone(nan_UF), clone(nan_UF));

    if (isInfinity_UF(x->Re) && isZero_UF(x->Im))
        return clone_Complex(x);

    if (isInfinity_UF(x->Re) && isNormal_UF(x->Im))
        return create_UFComplex(mul_UF(infinity_UF, cos_UF(x->Im)),
                                      mul_UF(infinity_UF, sin_UF(x->Im)));
    
    if (isNan_UF(x->Re) && isZero_UF(x->Im))
        return clone_Complex(x);

    if (isNan_UF(x->Re) && !isNan_UF(x->Im) && !isZero_UF(x->Im))
        return create_UFComplex(clone(nan_UF), clone(nan_UF));

    if (isNan_UF(x->Re) && isNan_UF(x->Im))
        return create_UFComplex(clone(nan_UF), clone(nan_UF));

    real = exp_UF(x->Re);
    real = sub_UF(real, div_UF(One, real));
    real = mul_UF(real, cos_UF(x->Im));
    real->exp--;

    imag = exp_UF(x->Re);
    imag = add_UF(imag, div_UF(One, imag));
    imag = mul_UF(imag, sin_UF(x->Im));
    imag->exp--;

    res = create_UFComplex(real, imag);

    return res;
}

UnifloatComplex* ctanh_UF(UnifloatComplex* x)
{
    UnifloatComplex* res;
    Unifloat* Zero = createZero_UF();
    Unifloat* Two = convertInteger_UF(2);

    if (x->Re->sign == -1)
    {
        x->Re->sign = 1;
        x->Im->sign = (-1) * x->Im->sign;
        res = ctanh_UF(x);
        x->Re->sign = -1;
        x->Im->sign = (-1) * x->Im->sign;
        res->Im->sign = (-1) * res->Im->sign;
        res->Re->sign = (-1) * res->Re->sign;
        return res;
    }

    if (x->Im->sign == -1)
    {
        x->Im->sign = 1;
        res = ctanh_UF(x);
        x->Im->sign = -1;
        res->Im->sign = (-1) * res->Im->sign;
        return res;
    }    

    if (isZero_UF(x->Re) && isZero_UF(x->Im))
        return clone_Complex(x);

    if (isNormal_UF(x->Re) && isInfinity_UF(x->Im))
        return create_UFComplex(clone(nan_UF), clone(nan_UF));

    if (isNormal_UF(x->Re) && isNan_UF(x->Im))
        return create_UFComplex(clone(nan_UF), clone(nan_UF));

    if (isInfinity_UF(x->Re) && isNormal_UF(x->Im))
        return create_UFComplex(createOne_UF(), mul_UF(Zero, sin_UF(mul_UF(x->Im, Two))));

    if (isInfinity_UF(x->Re) && isInfinity_UF(x->Im))
        return create_UFComplex(createOne_UF(), createZero_UF());

    if (isInfinity_UF(x->Re) && isNan_UF(x->Im))
        return create_UFComplex(createOne_UF(), createZero_UF());

    if (isNan_UF(x->Re) && isZero_UF(x->Im))
        return create_UFComplex(clone(nan_UF), createZero_UF());

    if (isNan_UF(x->Re) && !isNan_UF(x->Im) && !isZero_UF(x->Im))
        return create_UFComplex(clone(nan_UF), clone(nan_UF));

    if (isNan_UF(x->Re) && isNan_UF(x->Im))
        return create_UFComplex(clone(nan_UF), clone(nan_UF));

    return div_UFComplex(csinh_UF(x), ccosh_UF(x));   
}
