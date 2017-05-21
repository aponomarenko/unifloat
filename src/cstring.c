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

#include <string.h>
#include <stdlib.h>
#include <stdio.h> //perror

#include "unifloat/cstring.h"

#ifdef DEBUG_ON
long MEM = 0;
#endif

CString* create_CString(const char* buf)
{
    int len = strlen(buf);
    CString* cstr = (CString*) malloc(sizeof(CString));         
    if(cstr == NULL)
    {
        perror("Malloc error on creating a CString");
        exit(-1);
    }

    #ifdef DEBUG_ON
    MEM += sizeof(CString);
    #endif

    cstr->len = len;
    cstr->buf = (char*) malloc(len * sizeof(char));

    if(cstr->buf == NULL)
    {
        perror("Malloc error on creating a CString data");
        exit(-1);
    }

    #ifdef DEBUG_ON    
    MEM += len * sizeof(char);
    #endif
    memcpy(cstr->buf, buf, len);
    
    return cstr;
}

void delete_CString(CString* cstr)
{
    #ifdef DEBUG_ON    
    MEM -= cstr->len;
    #endif

    free(cstr->buf);

    #ifdef DEBUG_ON    
    MEM -= sizeof(CString);
    #endif

    free(cstr);
}

CString* clone_CString(CString* cstr)
{
    CString* new_cstr = (CString*) malloc(sizeof(CString)); 
    if(new_cstr == NULL)
    {
        perror("Malloc error on cloning a CString");
        exit(-1);
    }

    #ifdef DEBUG_ON
    MEM += sizeof(CString);
    #endif

    new_cstr->len = cstr->len;
    new_cstr->buf = (char*) malloc(cstr->len * sizeof(char));
    if(new_cstr->buf == NULL)
    {
        perror("Malloc error on cloning a CString data");
        exit(-1);
    }
    
    #ifdef DEBUG_ON
    MEM += cstr->len * sizeof(char);
    #endif

    memcpy(new_cstr->buf, cstr->buf, cstr->len);

    return new_cstr;
}

int length_CString(CString* cstr) {
    return cstr->len;
}

char charAt_CString(CString *cstr, int i) {
    return cstr->buf[i];
}

CString* substring_CString(CString* src, int start, int len)
{
    CString* cstr = (CString*) malloc(sizeof(CString));

    #ifdef DEBUG_ON
    MEM += sizeof(CString);
    #endif

    if(start + len > length_CString(src))
        len = length_CString(src) - start;

    cstr->len = len;
    cstr->buf = (char*) malloc(len * sizeof(char));

    #ifdef DEBUG_ON
    MEM += len * sizeof(char);
    #endif

    memcpy(cstr->buf, &src->buf[start], len);
    return cstr;
}

int indexOfChar_CString(CString* cstr, char c)
{
    int i;
    for(i=0; i<cstr->len; i++)
    {
        if(cstr->buf[i] == c)
            return i;
    }
    return -1;
}
