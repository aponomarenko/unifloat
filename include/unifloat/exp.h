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

/**\file exp.h
*\brief Real exponential functions. */

#ifndef EXP_H_
#define EXP_H_

#include "unifloat/unifloat.h"

Unifloat* exp_UF(Unifloat* x);
Unifloat* exp2_UF(Unifloat* x);
Unifloat* expm1_UF(Unifloat* x);
Unifloat* log_UF(Unifloat* x);
Unifloat* log1p_UF(Unifloat* x);
Unifloat* sqrt_UF(Unifloat* x);

#endif //EXP_H_
