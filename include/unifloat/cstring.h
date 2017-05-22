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

/**\file cstring.h
*\brief Auxiliary structure and functions to work with text. */

#ifndef CSTRING_H_
#define CSTRING_H_

#include "unifloat/config.h"

/*!
\brief The structure that represents a string in libunifloat.
*/
typedef struct CString 
{
    int len;
    char* buf;
} CString;

#ifdef DEBUG_ON
extern long MEM;
#endif

void delete_CString(CString* cstr);

CString* create_CString(const char* buf);
CString* clone_CString(CString* cstr);
int length_CString( CString* cstr);
char charAt_CString(CString *cstr, int i);
CString* substring_CString(CString* src, int start, int len);

int indexOfChar_CString (CString* cstr, char c);

#endif //CSTRING_H_
