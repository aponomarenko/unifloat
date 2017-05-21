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

/**\file bessel.h
*\brief Bessel functions.*/

#ifndef BESSEL_H_
#define BESSEL_H_

#include "unifloat/unifloat.h"

Unifloat* j0_UF(Unifloat* x);
Unifloat* j1_UF(Unifloat* x);
Unifloat* jn_UF(int n, Unifloat* x);

Unifloat* jnPowerSeries_UF(int n, Unifloat* x);
Unifloat* jnHankel_UF(int n, Unifloat* x);
Unifloat* jnSteed_UF(int n, Unifloat* x);
Unifloat* jnMeisselFirst_UF(int n, Unifloat* x);
Unifloat* jnMeisselSecond_UF(int n, Unifloat* x);
Unifloat* jnRecurrent_UF(int n, Unifloat* x);

Unifloat* y0_UF(Unifloat* x);
Unifloat* y1_UF(Unifloat* x);
Unifloat* yn_UF(int n, Unifloat* x);

Unifloat* ynPowerSeries_UF(int n, Unifloat* x);
Unifloat* ynHankel_UF(int n, Unifloat* x);
Unifloat* ynSteed_UF(int n, Unifloat* x);
Unifloat* ynMeisselFirst_UF(int n, Unifloat* x);
Unifloat* ynMeisselSecond_UF(int n, Unifloat* x);
Unifloat* ynRecurrent_UF(int n, Unifloat* x);

int besselMethod_UF(int n, Unifloat* x);

#endif// BESSEL_H_
