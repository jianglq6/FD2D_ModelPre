/***************************************************************************
 *
 * This function is used for the source term
 *
 * Authors: Luqian Jiang   Email: jlq0608@gmail.com
 *          Wei Zhang      Email: zhangwei@sustc.edu.cn
 *
 * Copyright (c) 2018 - 2019 zwlab
 *
 * History: 08/2019: Original version created by Luqian Jiang
 *
 ***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


// gauss and it deriv.
float fun_gauss(float t, float a, float t0)
{
    float f;
    if(fabs(t0)>SEIS_ZERO && (t<=0.0 || t>=2*t0)) {
        f = 0;
    } else {
        f = exp(-(t-t0)*(t-t0)/(a*a))/(sqrt(PI)*a);
    }

    return f;
}

float fun_gauss_deriv(float t, float a, float t0)
{
    float f;
    if(fabs(t0)>SEIS_ZERO && (t<=0.0 || t>=2*t0)) {
        f = 0;
    } else {
        f = exp(-(t-t0)*(t-t0)/(a*a))/(sqrt(PI)*a)*(-2*(t-t0)/(a*a));
    }

    return f;
}

// ricker and it deriv.
float fun_ricker(float t, float fc, float t0)
{
    float u, f0,v;
    if(t<=0.0) {
        v = 0.0;
    } else {
        f0 = sqrt(PI)/2.0;
        u = (t-t0)*2.0*PI*fc;
        v = (u*u/4-0.5)*exp(-u*u/4)*f0;
    }

    return v;
}

float fun_ricker_deriv(float t, float fc, float t0)
{
    float u, f0,v;
    if(t<=0.0) {
        v = 0.0;
    } else {
        f0 = sqrt(PI)/2.0;
        u = (t-t0)*2.0*PI*fc;
        v = u*(1.5-u*u/4)*exp(-u*u/4)*f0*PI*fc;
    }

    return v;
}

// bell, deriv and it integer.
float fun_bell(float t, float riset)
{
    float v;
    if(t>0.0 && t<riset) {
        v = (float)((1.0-cos(2*PI*t/riset))/riset);
    } else {
        v = 0;
    }
    return v;
}

float fun_bell_deriv(float t, float riset)
{
    float v;
    if(t>0.0 && t<riset) {
        v = (float)(2*PI*sin(2*PI*t/riset)/riset);
    } else {
        v = 0;
    }
    return v;
}

float fun_bell_int(float t, float riset)
{
    float v;
    if(t<=0.0) {
        v = 0;
    } else if(t<riset) {
        v = (float)(t/riset - sin(2*PI*t/riset)/(2*PI));
    } else {
        v = 1;
    }

    return v;
}

float fun_triangle(float t, float riset)
{
    float t0, v;
    t0 = riset/2;

    if(t>riset) {
        v = 0;
    } else if(t>t0) {
        v = 2/t0-t/(t0*t0);
    } else if(t>0) {
        v = t/(t0*t0);
    } else {
        v=0;
    }

    return v;
}

float fun_triangle_int(float t, float riset)
{
    float t0, v;
    t0 = riset/2.0;

    if(t>riset) {
        v = 1.0;
    } else if(t>t0) {
        v = -0.5*t*t/(t0*t0) + 2*t/t0 - 1;
    } else if(t>0) {
        v = 0.5*t*t/(t0*t0);
    } else {
        v = 0.0;
    }

    return v;
}

float fun_bshift(float t, float riset, float t0)
{
    float bshift, v;

    bshift = riset/2;
    if(t>0 && t<riset) {
        v = 0.5/t0/(cosh((t-bshift)/t0)*(t-bshift)/t0);
    } else {
        v = 0;
    }

    return v;
}

float fun_delta(float t, float t0)
{
    float v;

    if(t>=0 && t<t0) {
        v = 1.0/t0;
    } else {
        v = 0;
    }

    return v;
}



