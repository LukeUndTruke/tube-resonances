#ifndef GEOMETRIES_H
#define GEOMETRIES_H

// inline constexpr double c = 343.;
// inline constexpr double rho0 = 1.225;


#include<vector>
#include<complex>
#include "medium.h"
#include "math_helpers.h"

struct Geometry {
    virtual std::complex<double> inputImpedance(double freq, std::complex<double> Z_end) const = 0;
    struct TransferMatrix{
        std::complex<double> t11, t12, t21, t22;
    };
    virtual TransferMatrix T(double x, double freq) const = 0;
    virtual void computeFieldFromMatrix(double freq,
        std::complex<double> Z_end, 
        std::vector<double>& x, 
        std::vector<std::complex<double>>& p, 
        std::vector<std::complex<double>>&u) const = 0;
    
    virtual ~Geometry(){};
};

struct Cylinder : public Geometry {
    double L;  //length
    double R;  //radius
    double Z0;

    Cylinder(double L_, double R_): L(L_), R(R_), Z0(acoustics::rho0*acoustics::c/(constants::PI*R*R)) {};

    std::complex<double> inputImpedance(double freq, std::complex<double> Z_end) const override;

    //optional explicity analytic solution for cylinder
    void computeAnalyticField(double freq, 
        std::complex<double> Z_end, 
        std::vector<double>& x, 
        std::vector<std::complex<double>>& p, 
        std::vector<std::complex<double>>& u) const;

    TransferMatrix T(double x, double freq) const override;

    void computeFieldFromMatrix(double freq,
        std::complex<double> Z_end, 
        std::vector<double>& x, 
        std::vector<std::complex<double>>& p, 
        std::vector<std::complex<double>>&u) const override;
};

struct Cone : public Geometry {
    double L;  //slant length
    double R1; //left radius
    double R2; //right radius

    Cone(double L_, double R1_, double R2_): L(L_), R1(R1_), R2(R2_) {};

    std::complex<double> inputImpedance(double freq, std::complex<double> Z_end) const override;

    static TransferMatrix coneMatrix(double length, double r_start, double r_end, double freq);
    TransferMatrix T(double x, double freq) const override;

    void computeFieldFromMatrix(double freq,
        std::complex<double> Z_end, 
        std::vector<double>& x, 
        std::vector<std::complex<double>>& p, 
        std::vector<std::complex<double>>&u) const override;
};

#endif 