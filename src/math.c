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

#include "unifloat/unifloat.h"
#include "unifloat/trig.h"
#include "unifloat/exp.h"
#include "unifloat/bessel.h"
#include "unifloat/gamma.h"

static double call_UF_I(caller_UF func, double x)
{
    double res;
    
    initialize_UF();
    res = call_UF(func, x);
    finalize_UF();
    
    return res;
}

static double call_UF_NX_I(caller_UF_NX func, int n, double x)
{
    double res;
    
    initialize_UF();
    res = call_UF_NX(func, n, x);
    finalize_UF();
    
    return res;
}

//SIN
long double sinl(long double x) {
    return call_UF_I(&sin_UF, x);
}

double sin(double x) {
    return call_UF_I(&sin_UF, x);
}

float sinf(float x) {
    return call_UF_I(&sin_UF, x);
}

//ASIN
long double asinl(long double x) {
    return call_UF_I(&asin_UF, x);
}

double asin(double x) {
    return call_UF_I(&asin_UF, x);
}

float asinf(float x) {
    return call_UF_I(&asin_UF, x);
}

//COS
long double cosl(long double x) {
    return call_UF_I(&cos_UF, x);
}

double cos(double x) {
    return call_UF_I(&cos_UF, x);
}

float cosf(float x) {
    return call_UF_I(&cos_UF, x);
}

//ACOS
long double acosl(long double x) {
    return call_UF_I(&acos_UF, x);
}

double acos(double x) {
    return call_UF_I(&acos_UF, x);
}

float acosf(float x) {
    return call_UF_I(&acos_UF, x);
}

//TAN
long double tanl(long double x) {
    return call_UF_I(&tan_UF, x);
}

double tan(double x) {
    return call_UF_I(&tan_UF, x);
}

float tanf(float x) {
    return call_UF_I(&tan_UF, x);
}

//ATAN
long double atanl(long double x) {
    return call_UF_I(&atan_UF, x);
}

double atan(double x) {
    return call_UF_I(&atan_UF, x);
}

float atanf(float x) {
    return call_UF_I(&atan_UF, x);
}

//EXP
long double expl(long double x) {
    return call_UF_I(&exp_UF, x);
}

double exp(double x) {
    return call_UF_I(&exp_UF, x);
}

float expf(float x) {
    return call_UF_I(&exp_UF, x);
}

//LOG
long double logl(long double x) {
    return call_UF_I(&log_UF, x);
}

double log(double x) {
    return call_UF_I(&log_UF, x);
}

float logf(float x) {
    return call_UF_I(&log_UF, x);
}

//J0
long double j0l(long double x) {
    return call_UF_I(&j0_UF, x);
}

double j0(double x) {
    return call_UF_I(&j0_UF, x);
}

float j0f(float x) {
    return call_UF_I(&j0_UF, x);
}

//J1
long double j1l(long double x) {
    return call_UF_I(&j1_UF, x);
}

double j1(double x) {
    return call_UF_I(&j1_UF, x);
}

float j1f(float x) {
    return call_UF_I(&j1_UF, x);
}

//Jn
long double jnl(int n, long double x) {
    return call_UF_NX_I(&jn_UF, n, x);
}

double jn(int n, double x) {
    return call_UF_NX_I(&jn_UF, n, x);
}

float jnf(int n, float x) {
    return call_UF_NX_I(&jn_UF, n, x);
}

//Y0
long double y0l(long double x) {
    return call_UF_I(&y0_UF, x);
}

double y0(double x) {
    return call_UF_I(&y0_UF, x);
}

float y0f(float x) {
    return call_UF_I(&y0_UF, x);
}

//Y1
long double y1l(long double x) {
    return call_UF_I(&y1_UF, x);
}

double y1(double x) {
    return call_UF_I(&y1_UF, x);
}

float y1f(float x) {
    return call_UF_I(&y1_UF, x);
}

//Yn
long double ynl(int n, long double x) {
    return call_UF_NX_I(&yn_UF, n, x);
}

double yn(int n, double x) {
    return call_UF_NX_I(&yn_UF, n, x);
}

float ynf(int n, float x) {
    return call_UF_NX_I(&yn_UF, n, x);
}

//GAMMA
long double gammal(long double x) {
    return call_UF_I(&gamma_UF, x);
}

double gamma(double x) {
    return call_UF_I(&gamma_UF, x);
}

float gammaf(float x) {
    return call_UF_I(&gamma_UF, x);
}

//TGAMMA
long double tgammal(long double x) {
    return call_UF_I(&tgamma_UF, x);
}

double tgamma(double x) {
    return call_UF_I(&tgamma_UF, x);
}

float tgammaf(float x) {
    return call_UF_I(&tgamma_UF, x);
}

//LGAMMA
long double lgammal(long double x) {
    return call_UF_I(&lgamma_UF, x);
}

double lgamma(double x) {
    return call_UF_I(&lgamma_UF, x);
}

float lgammaf(float x) {
    return call_UF_I(&lgamma_UF, x);
}
