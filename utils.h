#ifndef UTILS_H_INCLUDED
#define UTILS_H_INCLUDED

#include <iostream>
#include <vector>
#include <cmath>
#include <numbers>
#include <complex>
#include "polynomial.h"
#include "config.h"

using val_t = long double;

const val_t PI    = std::numbers::pi_v<val_t>;
const val_t EPS   = 1e-17;


std::vector<val_t> find_roots(const Polynomial<val_t>&, val_t, val_t);
val_t root_binary_search(const Polynomial<val_t>& p, val_t l, val_t r);

#endif 