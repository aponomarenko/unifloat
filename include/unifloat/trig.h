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

/**\file trig.h
*\brief Real trigonometric functions. */

#ifndef TRIG_H_
#define TRIG_H_

#include "unifloat/unifloat.h"

/**\brief The sine of the given Unifloat number.
*
*\return The sine of \b x.
*\n The function creates and returns the object that have to be removed later using \b delete_UF. */
Unifloat* sin_UF(Unifloat* x);

/**\brief The cosine of the given Unifloat number.
*
*\return The cosine of \b x.
*\n The function creates and returns the object that have to be removed later using \b delete_UF. */
Unifloat* cos_UF(Unifloat* x);

/**\brief The tangent of the given Unifloat number.
*
*\return The tan of \b x.
*\n The function creates and returns the object that have to be removed later using \b delete_UF. */
Unifloat* tan_UF(Unifloat* x);

/**\brief The arc sine of the given Unifloat number.
*
*\return The arc sine of \b x.
*\n The function creates and returns the object that have to be removed later using \b delete_UF. */
Unifloat* asin_UF(Unifloat* x);

/**\brief The arc cosine of the given Unifloat number.
*
*\return The arc cosine of \b x.
*\n The function creates and returns the object that have to be removed later using \b delete_UF. */
Unifloat* acos_UF(Unifloat* x);

/**\brief The arc tangent of the given Unifloat number.
*
*\return The arc tangent of \b x.
*\n The function creates and returns the object that have to be removed later using \b delete_UF. */
Unifloat* atan_UF(Unifloat* x);

/**\brief The arc tangent of Unifloat number divided by another one.
*
*\return The arc tangent of <b>x/y</b> .
*\n The function creates and returns the object that have to be removed later using \b delete_UF. */
Unifloat* atan2_UF(Unifloat* x, Unifloat* y);

/**\brief Make a number to be between 0 and 2*Pi.
*
*\return The function creates and returns the object that have to be removed later using \b delete_UF. */
Unifloat* arrangeArgument_UF(Unifloat* x);

/**\brief Auxiliary function needed to calculate atan iteratively.
*
*\return The function creates and returns the object that have to be removed later using \b delete_UF. */
Unifloat* calcAtan_UF(Unifloat* x, Unifloat* y, int i);

/**\brief Auxiliary function needed to calculate tan iteratively.
*
*\return The function creates and returns the object that have to be removed later using \b delete_UF. */
Unifloat* calcTan_UF(Unifloat* x, int i);

#endif //TRIG_H_
