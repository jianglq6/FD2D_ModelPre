#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PI 3.14159265

/*  the zero-order modified Bessel function of the first kind */
double first_modified_Bessel(int n,double x)
{
    int i,m;
    double t,y,p,b0,b1,q;
    static double a[7]={ 1.0,3.5156229,3.0899424,1.2067492,
                0.2659732,0.0360768,0.0045813};
    static double b[7]={ 0.5,0.87890594,0.51498869,
                0.15084934,0.02658773,0.00301532,0.00032411};
    static double c[9]={ 0.39894228,0.01328592,0.00225319,
                -0.00157565,0.00916281,-0.02057706,
                0.02635537,-0.01647633,0.00392377};
    static double d[9]={ 0.39894228,-0.03988024,-0.00362018,
                 0.00163801,-0.01031555,0.02282967,
                -0.02895312,0.01787654,-0.00420059};
    if ( n < 0 )
        n = -n;
    t = fabs(x);
    if ( n != 1 )
    {
        if ( t < 3.75 )
        {
            y = (x/3.75)*(x/3.75);
            p = a[6];
            for (i = 5; i >= 0; i--)
                p = p*y + a[i];
        }
        else
        {
            y = 3.75/t;
            p = c[8];
            for (i = 7; i >= 0; i--)
                p = p*y+c[i];
            p = p*exp(t)/sqrt(t);
        }
    }
    if ( n == 0 )
        return(p);
    q = p;
    if ( t < 3.75 )
    {
        y = (x/3.75)*(x/3.75);
        p = b[6];
        for (i=5; i>=0; i--)
            p = p*y+b[i];
        p = p*t;
    }
    else
    {
        y = 3.75/t;
        p = d[8];
        for (i=7; i>=0; i--)
            p = p*y+d[i];
        p = p*exp(t)/sqrt(t);
    }
    if (x < 0.0)
        p = -p;
    if (n == 1)
        return(p);
    if (x == 0.0)
        return(0.0);
    y = 2.0/t;
    t = 0.0; b1 = 1.0; b0 = 0.0;
    m = n+(int)sqrt(40.0*n);
    m = 2*m;

    for (i = m; i > 0; i--)
    {
        p = b0+i*y*b1; b0=b1; b1=p;
        if (fabs(b1)>1.0e+10)
        {
            t  = t*1.0e-10;
            b0 = b0*1.0e-10;
            b1 = b1*1.0e-10;
        }
        if (i == n)
            t = b0;
    }
    p = t*q/b1;
    if ((x<0.0)&&(n%2==1))
        p = -p;
    return(p);
}
