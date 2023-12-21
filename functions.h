#ifndef FUNCTIONS_H_INCLUDED
#define FUNCTIONS_H_INCLUDED

#include <functional>
#include <algorithm>
#include <stdio.h>
#include <cmath>
#include <ctime>
#include "matrix.h"
#include "polynomial.h"
#include "utils.h"

val_t fp(val_t);
val_t f(val_t);
val_t p(val_t);


val_t integrate_qaws(const std::function<val_t(val_t)>& f, val_t a, val_t b, val_t alpha, val_t beta, val_t h, bool gauss = 1, size_t k = 3);

val_t integrate_newton(const std::function<val_t(val_t)>& f, val_t a, val_t b, val_t beta, size_t n);

val_t integrate_gauss(const std::function<val_t(val_t)>& f, val_t a, val_t b, val_t beta, size_t n);

#endif

