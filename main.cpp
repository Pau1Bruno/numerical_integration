#include "functions.h"
#include <string>
#include <fstream>
#include <ios>


signed main()
{
    using namespace std;
    std::cout.precision(16);
    std::cout.setf(std::ios::fixed);
   
    val_t a = val_t(7)/10, b = val_t(32)/10, beta = val_t(1)/4;
    bool gauss = true;

    int n = 2;
    val_t h = (b - a);
    val_t m,R12,R21,R23,R32;
    val_t L = 1.5;
    val_t sum1, sum2 = integrate_qaws(f,a,b,0,beta,h,gauss,n);
    h /= L;
    // val_t sum2 = integrate_qaws(f,a,b,0,beta,h,gauss,n);
    // h /= L;
    val_t sum3 = integrate_qaws(f,a,b,0,beta,h,gauss,n);

    // m = abs(log(abs(sum3-sum2)/abs(sum2-sum1)) / log(L));
    // R12 = abs(sum2-sum1)/(1 - pow(L,-m));
    // R21 = abs(sum2-sum1)/(pow(L,m)-1);
    // R23 = abs(sum3-sum2)/(1 - pow(L,-m));
    // R32 = abs(sum3-sum2)/(pow(L,m)-1);

    // cout << "ИКФ " << (gauss ? "Гаусса " : "Ньютона-Котса ") << sum3 << "\n";
    // cout << "R12 = " << R12 << "\n";
    // cout << "R21 = " << R21 << "\n";
    // cout << "R23 = " << R23 << "\n";
    // cout << "R32 = " << R32 << "\n";

    for (int i = 0; ; i++)
    {
        h /= L;
        sum1 = sum2;
        sum2 = sum3;
        sum3 = integrate_qaws(f,a,b,0,beta,h,gauss,n);
        m = abs(log(abs(sum3-sum2)/abs(sum2-sum1)) / log(L));
        
        R12 = abs(sum2-sum1)/(1 - pow(L,-m));
        R21 = abs(sum2-sum1)/(pow(L,m)-1);
        R23 = abs(sum3-sum2)/(1 - pow(L,-m));
        R32 = abs(sum3-sum2)/(pow(L,m)-1);
        
        cout << "h = " << h << "\n" << " m = " << m << "\n";
        cout << "ИКФ " << (gauss ? "Гаусса " : "Ньютона-Котса ") << sum3 << "\n";
        cout << "R12 = " << R12 << "\n";
        cout << "R21 = " << R21 << "\n";
        cout << "R23 = " << R23 << "\n";
        cout << "R32 = " << R32 << "\n";
    }
    
    return 0;
}