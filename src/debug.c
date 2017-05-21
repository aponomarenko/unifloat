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

#include <stdio.h>

#include "unifloat/config.h"

#ifdef DEBUG_ON

#include "unifloat/unifloat.h"
#include "unifloat/debug.h"

static long inuse_start_UF = 0;

void PRINT_InUseComment(const char* comment) {
    printf("Unifloats: %ld (%s)\n", count_UF, comment);
}

void PRINT_InUse() {
    printf("Unifloats: %ld\n", count_UF);
}

void PRINT_InUseStart()
{
    printf("Unifloats at start: %ld\n", count_UF);
    inuse_start_UF = count_UF;
}

void PRINT_InUseEnd() {
    printf("Unifloats at end: %ld, delta=%ld\n", count_UF, count_UF - inuse_start_UF);
}

void PRINT_U(const char* str, Unifloat* u)
{
    double d = convertUnifloat_Double(u);
    printf("%s: %f\n", str, d);
}

void PRINT_C(const char* str) {
    printf("%s\n", str);
}

void PRINT_I(unsigned int u)
{
    int i, j;
    printf("I ");
    for(i=3; i>=0; i--)
    {
        for(j=7; j>=0; j--)
            printf("%d", (u>>(j + i*8)) % 2);
        printf(" ");
    }
    printf("\n");
}

void PRINT_ALL(const char* str, Unifloat* u)
{
    //PRINT_U(str,u);
    printf("%s:\n", str);
    PRINT_I(u->mant[0]);
    PRINT_I(u->mant[1]);
    PRINT_I(u->mant[2]);
    PRINT_I(u->exp);
}

#endif //DEBUG_ON
