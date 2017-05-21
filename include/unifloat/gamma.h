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

/**\file gamma.h
*\brief Gamma functions.*/

#ifndef GAMMA_H_
#define GAMMA_H_

#include "unifloat/unifloat.h"

int signgam;
Unifloat* gamma_UF(Unifloat* x);
Unifloat* lgamma_UF(Unifloat* x);
Unifloat* tgamma_UF(Unifloat* x);
Unifloat* gammaSeries_UF(Unifloat* x);

#endif //GAMMA_H_
