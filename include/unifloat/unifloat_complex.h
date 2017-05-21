/*
    This file contains definition of the UnifloatComplex data structure
    and basic operations.
    
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

/**\file unifloat_complex.h
*\brief Basic UnifloatComplex operations.*/

#ifndef UNIFLOAT_COMPLEX_H_
#define UNIFLOAT_COMPLEX_H_

#include "unifloat/unifloat.h"

/**\brief Represents a complex number in libunifloat.
*
*The main structure for complex calculations.*/
typedef struct UnifloatComplex
{
    /**\brief The real part of complex number*/
    Unifloat* Re;
    /**\brief The imaginary part of complex number*/
    Unifloat* Im;
} UnifloatComplex;

/**\brief Creates a copy of given UnifloatComplex number.
*
*\return The copy of given UnifloatComplex number.
*\n The function creates and returns the object that have to be removed later using \b delete_UFComplex.*/
UnifloatComplex* clone_Complex(UnifloatComplex* src);

/**\brief Copies a given Unifloat object data to another one.*/
void copy(Unifloat* src, Unifloat* dst);

/**\brief Copies a given UnifloatComplex object data to another one.*/
void copy_Complex(UnifloatComplex* src, UnifloatComplex* dst);

/**\brief Create a new UnifloatComplex number.
*
*Mantissa is filled with zeros. 
*\return The new UnifloatComplex object. 
*\n The created object have to be removed later using \b delete_UFComplex.*/
UnifloatComplex* create_UFComplex( Unifloat* Re, Unifloat* Im);

/**\brief Deletes the UnifloatComplex number.
*
*This function have to be called for every created UnifloatComplex to free the memory it uses.*/
void delete_UFComplex( UnifloatComplex* x);

/**\brief Deletes the list of UnifloatComplex numbers.*/
void delete_UFsComplex(UnifloatComplex* p1, ...);


/**\brief Get the absolute value of the given UnifloatComplex number.
*
*\return The absolute value of \b x.
*\n The function creates and returns the object that have to be removed later using \b delete_UFComplex.*/
Unifloat* abs_UFComplex(UnifloatComplex* x);

/**\brief Get the argument of the given UnifloatComplex number.
*
*\return The argument of \b x.
*\n The function creates and returns the object that have to be removed later using \b delete_UFComplex.*/
Unifloat* carg_UF(UnifloatComplex* x);

/**\brief The sum of two Unifloat numbers.
*
*\return The UnifloatComplex object that represents a sum of \b x and \b y.
*\n The function creates and returns the object that have to be removed later using \b delete_UFComplex. */
UnifloatComplex* add_UFComplex(UnifloatComplex* x, UnifloatComplex* y);

/**\brief The difference between two Unifloat numbers.
*
*\return The UnifloatComplex object that represents a difference between \b x and \b y.
*\n The function creates and returns the object that have to be removed later using \b delete_UFComplex. */
UnifloatComplex* sub_UFComplex(UnifloatComplex* x, UnifloatComplex* y);

/**\brief The product of two Unifloat numbers.
*
*\return The UnifloatComplex object that represents a product of \b x and \b y.
*\n The function creates and returns the object that have to be removed later using \b delete_UFComplex. */
UnifloatComplex* mul_UFComplex(UnifloatComplex* x, UnifloatComplex* y);

/**\brief The ratio of two Unifloat numbers.
*
*\return The UnifloatComplex object that represents a ratio of \b x and \b y.
*\n The function creates and returns the object that have to be removed later using \b delete_UFComplex. */
UnifloatComplex* div_UFComplex(UnifloatComplex* x, UnifloatComplex* y);

#endif //UNIFLOAT_COMPLEX_H_
