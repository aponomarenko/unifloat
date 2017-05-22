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

/**\file debug.h
*\brief Auxiliary module to help programmer to debug the library. */

#ifndef UF_DEBUG_H_
#define UF_DEBUG_H_

#include "unifloat/config.h"

#ifdef DEBUG_ON
#include "unifloat/unifloat.h"

void PRINT_InUse(void);
void PRINT_InUseComment(const char* str);
void PRINT_InUseStart(void);
void PRINT_InUseEnd(void);
void PRINT_U(const char* str, Unifloat* u);
void PRINT_ALL(const char* str, Unifloat* u);
void PRINT_I(unsigned int I);
void PRINT_C(const char* str);

#define ST PRINT_InUseStart();
#define END PRINT_InUseEnd();

#endif //DEBUG_ON

#endif //UF_DEBUG_H_
