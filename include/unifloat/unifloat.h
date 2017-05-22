/*
    This file contains definition of the Unifloat data structure and
    basic operations.
    
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

/**\file unifloat.h
*\brief Basic Unifloat operations.*/

#ifndef UNIFLOAT_H_
#define UNIFLOAT_H_

#include "unifloat/config.h"
#include "unifloat/cstring.h"

/**\brief Maximal size of Unifloat mantissa*/
#define MAX_SIZE_UNIFLOAT 7
#define true 1
#define false 0
/**\brief Maximal count of mantissa bits, that can be used for calculations*/
#define PRECISION 90
/**\brief Specifies the accuracy of Unifloats comparing.*/
#define COMPARE_PRECISION 15

#define maxExp_DoubleT 1024
#define minExp_DoubleT -1021
#define digMant_DoubleT 53
#define digExp_DoubleT 11
#define size_DoubleT 8
#define sizeInLongs_DoubleT (sizeof(double) - 4) / sizeof(long) + 1

typedef unsigned int uint;
typedef unsigned long ulong;
typedef int Bool;

#ifdef DEBUG_ON
extern long count_UF;
#endif

/**Kind of Unifloat.*/
typedef enum UnifloatKind
{
    Infinity,
    NaN,
    Normal
} UnifloatKind;

/*!
\brief The main structure that represents a number in libunifloat.
*/
typedef struct Unifloat
{
    /**\brief The sign of Unifloat is the same as of this field.*/
    int  sign;
    /**\brief The power of two*/
    int  exp;
    /**\brief Mantissa of the number. 
    *
    *Not all bits are used. See \b PRECISION*/
    uint mant[MAX_SIZE_UNIFLOAT]; 
    /**\brief Kind of floating point number (Normal, Infinity, NaN).*/
    UnifloatKind kind;
} Unifloat;

/**\brief The number PI.*/
extern Unifloat* Pi;

/**\brief The number GAMMA.*/
extern Unifloat* Gamma;

/**\brief The number E.*/
extern Unifloat* E;

extern Unifloat* Ln_2;
extern Unifloat* Log2_e;
extern Unifloat* Log10_e;
extern Unifloat* Log2_10;
extern Unifloat* Epsilon;

/**\brief Initialize global variables.
*
*Call it when you start working with libunifloat.*/
void initialize_UF(void);

/**\brief Frees the memory used by global variables.
*
*Call it when you end working with libunifloat.*/
void finalize_UF(void);

/**\brief Used to test the Unifloat for overflow*/
extern Unifloat* max_UF;

/**\brief Used to test the Unifloat for underflow*/
extern Unifloat* min_UF;

/**\brief The positive infinity Unifloat*/
extern Unifloat* infinity_UF;

/**\brief The Not-A-Number Unifloat*/
extern Unifloat* nan_UF;

/**\brief Represents a template of any function that takes one Unifloat argument.*/
typedef Unifloat* (*caller_UF) (Unifloat*);

/**\brief Represents a template of any function that takes integer and Unifloat arguments.*/
typedef Unifloat* (*caller_UF_nx) (int, Unifloat*);

/**\brief Represents a template of any function that takes two Unifloat arguments.*/
typedef Unifloat* (*caller_UF_xy) (Unifloat*, Unifloat*);

/**\brief brief Function needed to simplify memory management.
*
* In code "x = Foo(x)" variable x shold be deleted first.
* Call "x = call1_arg1(&Foo, x)" and the old Unifloat x will be removed.*/
Unifloat* call1_arg1(caller_UF func, Unifloat* x);

/**\brief Function needed to simplify memory management.
*
* In code "x = Foo(x, y)" variable x shold be deleted first.
* Call "x = call2_arg1(&Foo, x, y)" and the old Unifloat x will be removed. */
Unifloat* call2_arg1(caller_UF_xy func, Unifloat* x, Unifloat* y);

/**\brief Function needed to simplify memory management.
*
* In code "x = Foo(y, x)" variable x shold be deleted first.
* Call "x = call2_arg2(&Foo, y, x)" and the old Unifloat x will be removed. */
Unifloat* call2_arg2(caller_UF_xy func, Unifloat* x, Unifloat* y);

/**\brief Create a new Unifloat number.
*
*Mantissa is filled with zeros. 
*\return The new Unifloat object. 
*\n The created object have to be removed later using \b delete_UF.*/
Unifloat* create_UF(int sign, int exponent, UnifloatKind kind);

/**\brief Delete a Unifloat number and frees the used memory.*/
void delete_UF(Unifloat* u);

/**\brief Delete a list of Unifloat numbers and frees the used memory.
*
*\n The last argument should be \b NULL*/
void delete_UFs(Unifloat* u, ...);

/**\brief Create a new Unifloat object that represent Zero.
*
*\return The new Unifloat Zero object.
*\n The created object have to be removed later using \b delete_UF.*/
Unifloat* createZero_UF(void);

/**\brief Create a new Unifloat object that represent One.
*
*\return The new Unifloat One object.
*\n The created object have to be removed later using \b delete_UF. */
Unifloat* createOne_UF(void);

/**\brief Creates a copy of given Unifloat number.
*
*\return The copy of given Unifloat number. 
*\n The function creates and returns the object that have to be removed later using \b delete_UF.*/
Unifloat* clone(Unifloat* src);

/**\brief Copies a given Unifloat number from /b src to /b dst.*/
void copy(Unifloat* src, Unifloat* dst);

/**\brief Check if the Kind of number is Normal.*/
Bool isNormal_UF(Unifloat* x);

/**\brief Check if the Kind of number is Infinity.*/
Bool isInfinity_UF(Unifloat* x);

/**\brief Check if the Kind of number is NaN.*/
Bool isNan_UF(Unifloat* x);

/**\brief Check if the Unifloat number is Zero.*/
Bool isZero_UF(Unifloat* x);

/**\brief Check if the Unifloat number is overflowed.*/
Bool isOverflow_UF(Unifloat* x);

/**\brief Check if the Unifloat number is underflowed.*/
Bool isUnderflow_UF(Unifloat* x);

/**\brief Returns a absolute value of given number.
*
*\return The Unifloat object that represents absolute value of a number \b x.
*\n Function creates and returns the object that have to be removed later using \b delete_UF.*/
Unifloat* abs_UF(Unifloat* x);

/**\brief Returns a normalized value of a given number.
*
*Normalized Unifloat have the first bit of mantissa set to 1.
*\return The Unifloat object that represents normalized value of a number \b x.
*\n The function creates and returns the object that have to be removed later using \b delete_UF.*/
Unifloat* normalize_UF(Unifloat* x);

/**\brief Sets one bit of Unifloat mantissa to a given value.
*
*\param index the bit number to be set.
*\param bit the value of bit*/
void setMant_UF(Unifloat* x, uint index, uint bit);

/**\brief Gets one bit of Unifloat mantissa.
*
*\param index the bit number to get.*/
uint getMant_UF(Unifloat* x, uint index);

/**\brief Rounds a Unifloat number. 
*
*\return The function creates and returns the object that have to be removed later using \b delete_UF.
*\param precision count of significant bits to save.*/
Unifloat* round_UF(Unifloat* x, uint precision);

/**\brief Rounds a Unifloat number. 
*
*\return The function creates and returns the object that have to be removed later using \b delete_UF.
*\param s if it is set to 0 - change the sign to opposite.
*Otherwise multiply the sign by s. */
Unifloat* changeSign_UF(Unifloat* x, int s);

/**\brief Compares two Unifloat numbers. 
*
*\return \b 0, if x==y (or they are both Infinity or NaN); \n \b 1, if x>y;
* \n \b -1, if x<y; \n \b 2, if one (and only one) of x or y is NaN.*/
int compare_UF(Unifloat* x, Unifloat* y);

/**\brief Checks that the given Unifloat number is Error.
*
*Error is a Unifloat with \b exp set to -150
*\return The result of comparing the internal Error object with the given one using \b compare_UF.*/
int compareWithError_UF(Unifloat* x);

/**\brief Compares two Unifloat numbers with the specified precision.
*
*\return \b 0, if x==y (or they are both Infinity or NaN); \n \b 1, if x>y;
* \n \b -1, if x<y; \n \b 2, if one (and only one) of x or y is NaN.*/
int compareWithPrecision_UF(Unifloat* x, Unifloat* y, int amount);

/**\brief The sum of two Unifloat numbers.
*
*\return The Unifloat object that represents a sum of \b x and \b y.
*\n The function creates and returns the object that have to be removed later using \b delete_UF. */
Unifloat* add_UF(Unifloat* x, Unifloat* y);

/**\brief The difference between two Unifloat numbers.
*
*\return The Unifloat object that represents a difference between \b x and \b y.
*\n The function creates and returns the object that have to be removed later using \b delete_UF. */
Unifloat* sub_UF(Unifloat* x, Unifloat* y);

/**\brief The product of two Unifloat numbers.
*
*\return The Unifloat object that represents a product of \b x and \b y.
*\n The function creates and returns the object that have to be removed later using \b delete_UF. */
Unifloat* mul_UF(Unifloat* x, Unifloat* y);

/**\brief The ratio of two Unifloat numbers.
*
*\return The Unifloat object that represents a ratio of \b x and \b y.
*\n The function creates and returns the object that have to be removed later using \b delete_UF. */
Unifloat* div_UF(Unifloat* x, Unifloat* y);

/**\brief Convert a given \b float number into Unifloat number.
*
*\return The Unifloat object that represents a given number \b x in Unifloat.
*\n The function creates and returns the object that have to be removed later using \b delete_UF. */
Unifloat* convertFloat_UF(float x);

/**\brief Convert a given \b double number into Unifloat number.
*
*\return The Unifloat object that represents a given number \b x in Unifloat.
*\n The function creates and returns the object that have to be removed later using \b delete_UF. */
Unifloat* convertDouble_UF(double x);

/**\brief Convert a given long double number into Unifloat number.
*
*\return The Unifloat object that represents a given number \b x in Unifloat.
*\n The function creates and returns the object that have to be removed later using \b delete_UF. */
Unifloat* convertLongDouble_UF(long double x);

/**\brief Convert a given Unifloat number into \b float number.
*
*\return \b float number that represents \b x. */
float convertUnifloat_Float(Unifloat* x);

/**\brief Convert a given Unifloat number into \b double number.
*
*\return \b double number that represents \b x. */
double convertUnifloat_Double(Unifloat* x);

/**\brief Convert a given Unifloat number into long double number.
*
*\return long double number that represents \b x. */
long double convertUnifloat_LongDouble(Unifloat* x);

/**\brief Convert a given \b integer number into Unifloat number.
*
*\return The Unifloat object that represents a given number \b number in Unifloat.
*\n The function creates and returns the object that have to be removed later using \b delete_UF. */
Unifloat* convertInteger_UF(int number);

/**\brief Convert a given Unifloat number into \b integer number.
*
*\param error The function returns an error code using this parameter. \n \b 0 means that it is possible to convert x to integer number,
*\n \b 1 means that conversion is impossible, \n \b 2 means that conversion isn't precisely and \n \b 3 means that result isn't representable in \b int type.
*\return The \b integer number that represent a given Unifloat number \b x. */
int convertUnifloat_Integer(Unifloat* x, int* error);

/**\brief Parse the \b CString that contains a binary number and create the Unifloat number.
*
*\return The Unifloat object that represents \b x in Unifloat.
*\n The function creates and returns the object that have to be removed later using \b delete_UF. */
Unifloat* convertBinaryString_UF(CString* number);

/**\brief Parse the \b CString that contains a decimal number and create the Unifloat number.
*
*\return The Unifloat object that represents \b x in Unifloat.
*\n The function creates and returns the object that have to be removed later using \b delete_UF. */
Unifloat* convertString_UF(CString* number);

/**\brief Raise the value of Unifloat number \b x to the power of \b n.
*
*\return The Unifloat object that represents \b x in the power of \b n in Unifloat.
*\n The function creates and returns the object that have to be removed later using \b delete_UF. */
Unifloat* power_UF(Unifloat* x, int n);

/**\brief Compute factorial of an integer number \b n.
*
*\return The Unifloat object that represents factorial of \b n in Unifloat.
*\n The function creates and returns the object that have to be removed later using \b delete_UF. */
Unifloat* factorial_UF(int n);

/**\brief Print value of the Unifloat number \b x in double format to stdout.*/
void print_UF(Unifloat* x);

/**\brief Call any Unifloat function \b func with the \b x argument in double.
*
*\return The return value of the function \b func in double. */
double call_UF(caller_UF func, double x);

/**\brief Call any Unifloat function \b func with the \b n and \b x arguments.
*
*\return The return value of the function \b func in double. */
double call_UF_nx(caller_UF_nx func, int n, double x);

#endif //UNIFLOAT_H_
