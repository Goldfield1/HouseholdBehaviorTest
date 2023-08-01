// functions related to utility and environment.
#ifndef MAIN
#define UTILS
#include "myheader.cpp"
#endif

namespace utils {
    double util(double c,double hours,par_struct *par) {
        double eta = par->eta;
        double beta = par->beta;
        double gamma = par->gamma;

        return pow(c,1+eta) / (1+eta) - beta * pow(hours,1+gamma) / (1+gamma);
    }

}