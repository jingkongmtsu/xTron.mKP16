/**
 * this file contains small utility functions for geom module 
 * author: fenglai liu
 */
#include<cmath>
#include "libgen.h"
#include "geomutil.h"
using namespace geomutil;

void geomutil::rotationByZ(Double& x, Double& y, Double angle) {
    Double tmpx;
    Double tmpy;
    Double cosgama = cos(angle);
    Double singama = sin(angle);

    // counter-clockwise
    tmpx = cosgama*x - singama*y;
    tmpy = singama*x + cosgama*y;
    x = tmpx;
    y = tmpy;
}

void geomutil::rotationByY(Double& x, Double& z, Double angle) {
    Double tmpx;
    Double tmpz;
    Double cosbeta = cos(angle);
    Double sinbeta = sin(angle);

    // counter-clockwise
    tmpx = cosbeta*x + sinbeta*z;
    tmpz =-sinbeta*x + cosbeta*z;
    x = tmpx;
    z = tmpz;
}

void geomutil::rotationByX(Double& y, Double& z, Double angle) {
    Double tmpy;
    Double tmpz;
    Double cosalfa = cos(angle);
    Double sinalfa = sin(angle);

    // counter-clockwise
    tmpy = cosalfa*y - sinalfa*z;
    tmpz = sinalfa*y + cosalfa*z;
    y = tmpy;
    z = tmpz;
}



