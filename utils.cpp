#include "utils.h"

std::vector<val_t> roots_quadratic(const Polynomial<val_t>& poly)
{
    val_t a = poly[2], b = poly[1], c = poly[0];
    val_t d = b*b - 4*a*c;

    if (d < 0)
        return {};

    return {(-b + sqrt(d))/2/a, (-b - sqrt(d))/2/a};
}

std::vector<val_t> find_roots(const Polynomial<val_t>& p, val_t l, val_t r)
{
    size_t n = p.deg();
    if (n == 0)
        return {};
    if (n == 1)
        return {-p[0]};
    if (n == 2)
        return roots_quadratic(p);

    return {};
} 