#include "functions.h"

val_t f (val_t x) 
{
    using namespace std;
    val_t res = 1.3 * cos(3.5 * x) * exp(2./3 * x) + 6 * sin(4.5 * x) * exp(-1./8 * x) + 5 * x;
    // cout << "f в точке " << x << " равна " << res << "\n";
    return res;
}

val_t moment(size_t k, val_t a, val_t b, val_t alpha)
{
    val_t power = k-alpha+1;
    return (pow(b, power) - pow(a, power)) / power;
}

std::vector<val_t> poly_interpol_quadrature(const std::vector<val_t>& nodes, val_t a, val_t b, val_t alpha)
{
    int n = nodes.size();

    Matrix<val_t> A = vandermonde(nodes, n-1, true);
    // std::cout << "Матрица: " << A << "\n";

    Row<val_t> c(n);
    for (int i = 0; i < n; i++)
        c[i] = moment(i,a,b,alpha);
    // std::cout << "Моменты: " << c << "\n";

    return LU_solve(A,c).getData();
}

std::vector<val_t> find_gaussian_nodes(val_t a, val_t b, val_t alpha, size_t n)
{
    std::vector<val_t> moments(2*n);
    for (size_t i = 0; i < 2*n; i++)
        moments[i] = moment(i,a,b,alpha);

    Matrix<val_t> M(n,n);
    for (size_t i = 0; i < n; i++)
        for (size_t j = 0; j < n; j++)
            M[i][j] = moments[i+j];

    std::vector<val_t> tmp = LU_solve(M,-Row(std::vector(moments.begin()+n, moments.begin()+2*n))).getData();

    tmp.push_back(1);
    std::reverse(tmp.begin(), tmp.end());
    Polynomial<val_t> node_poly = tmp;

    auto res = find_roots(node_poly, a, b);
    std::sort(res.begin(), res.end());

    return res;
}

val_t integrate_newton(const std::function<val_t(val_t)>& f, val_t a, val_t b, val_t alpha, size_t n)
{
    if (a > b)
        return 0;
    std::vector<val_t> nodes(n);
    val_t h = (b-a)/(n-1);
    for (size_t i = 0; i < n; i++)
            nodes[i] = a + h*i;

    
    std::vector<val_t> coeffs = poly_interpol_quadrature(nodes, a, b, alpha);
    for (val_t cf : coeffs)
        if (cf < 0)
        {
            std::cout << "Отрицательный вес\n";
            exit(-3);
        }

    val_t res = 0;
    for (size_t i = 0; i < n; i++)
            res += coeffs[i] * f(nodes[i]);

    return res;
}

val_t integrate_gauss(const std::function<val_t(val_t)>& f, val_t a, val_t b, val_t beta, size_t n)
{
    std::vector<val_t> nodes = find_gaussian_nodes(a, b, beta, n);
    
    if (nodes.size() != n)
    {
        std::cout << "Нет корней\n";
        exit(-1);
    }

    for (size_t i = 0; i < n; i++)
        if (nodes[i] < a || nodes[i] > b)
        {
            std::cout << "Корни находятся вне промежутка\n";
            exit(-2);
        }
        
    std::vector<val_t> coeffs = poly_interpol_quadrature(nodes, a, b, beta);

    for (val_t cf : coeffs)
        if (cf < 0)
        {
            std::cout << "Отрицательный вес\n";
            exit(-3);
        }

    val_t res = 0;
    for (size_t i = 0; i < n; i++)
            res += coeffs[i] * f(nodes[i]);

    return res;
}

static val_t integrate_qaws_wrapper(const std::function<val_t(val_t)>& f, val_t b, val_t beta, val_t h, bool gauss, size_t k)
{
    size_t n = ceil(double(b/h));
    h = b/n;

    val_t sum = 0, a = 0;
    for (size_t i = 0; i < n; i++)
    {
        sum += (gauss ? integrate_gauss(f, a, a+h, beta, k) : integrate_newton(f, a, a+h, beta, k));
        a += h;
    }

    return sum;
}

val_t integrate_qaws(const std::function<val_t(val_t)>& f, val_t a, val_t b, val_t alpha, val_t beta, val_t h, bool gauss, size_t k)
{
    if (alpha == 0)
        return integrate_qaws_wrapper([f,b](val_t x){return f(b-x);}, b-a, beta, h, gauss, k);

    if (beta == 0)
        return integrate_qaws_wrapper([f,a](val_t x){return f(x+a);}, b-a, alpha, h, gauss, k);

    return val_t(0)/0;
}