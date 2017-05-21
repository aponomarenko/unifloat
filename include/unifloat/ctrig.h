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

/**\file ctrig.h
*\brief Complex trigonometric functions. */

#ifndef CTRIG_H_
#define CTRIG_H_

#include "unifloat/unifloat_complex.h"
#include "unifloat/trig.h"

/**\brief The complex arc hyperbolic cosine of the given Unifloat number.
*
*\return The complex arc hyperbolic cosine of \b x.
*\n The function creates and returns the object that have to be removed later using \b delete_UF. */
UnifloatComplex* cacosh_UF(UnifloatComplex* x);

/**\brief The complex arc hyperbolic sine of the given Unifloat number.
*
*\return The complex arc hyperbolic sine of \b x.
*\n The function creates and returns the object that have to be removed later using \b delete_UF. */
UnifloatComplex* casinh_UF(UnifloatComplex* x);

/**\brief The complex arc hyperbolic sine of the given Unifloat number.
*
*\return The complex arc hyperbolic sine of \b x.
*\n The function creates and returns the object that have to be removed later using \b delete_UF. */
UnifloatComplex* catanh_UF(UnifloatComplex* x);

/**\brief The complex hyperbolic cosine of the given Unifloat number.
*
*\return The complex hyperbolic cosine of \b x.
*\n The function creates and returns the object that have to be removed later using \b delete_UF. */
UnifloatComplex* ccosh_UF(UnifloatComplex* x);

/**\brief The complex hyperbolic sine of the given Unifloat number.
*
*\return The complex hyperbolic sine of \b x.
*\n The function creates and returns the object that have to be removed later using \b delete_UF. */
UnifloatComplex* csinh_UF(UnifloatComplex* x);

/**\brief The complex hyperbolic tangent of the given Unifloat number.
*
*\return The complex hyperbolic tangent of \b x.
*\n The function creates and returns the object that have to be removed later using \b delete_UF. */
UnifloatComplex* ctanh_UF(UnifloatComplex* x);

#endif //CTRIG_H_
