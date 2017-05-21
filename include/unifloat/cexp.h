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

/**\file cexp.h
*\brief Complex exponential functions. */

#ifndef CEXP_H_
#define CEXP_H_

#include "unifloat/unifloat_complex.h"
#include "unifloat/exp.h"

UnifloatComplex* cexp_UF(UnifloatComplex* x);
UnifloatComplex* clog_UF(UnifloatComplex* x);
UnifloatComplex* clog10_UF(UnifloatComplex* x);
UnifloatComplex* cpow_UF(UnifloatComplex* x, UnifloatComplex* y);
UnifloatComplex* csqrt_UF(UnifloatComplex* x);

#endif //CEXP_H_
