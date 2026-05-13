#include "geometries.h"
#include <cmath>
#include <algorithm>


//Cylinder
std::complex<double> Cylinder::inputImpedance(double freq, std::complex<double> Z_end) const {
    // return 0;
    Geometry::TransferMatrix T_L = T(L, freq);
    return (T_L.t11*Z_end + T_L.t12) / (T_L.t21*Z_end + T_L.t22);
}

void Cylinder::computeAnalyticField(double freq, 
    std::complex<double> Z_end,
    std::vector<double>& x, 
    std::vector<std::complex<double>>& p, 
    std::vector<std::complex<double>>& u) 
    const {
        size_t N = x.size();
        if (p.size() != N) {
            p.resize(N);
        }
        if (u.size() != N) {
            u.resize(N);
        }
        const double k = 2.0*constants::PI*freq/acoustics::c;
        const std::complex<double> j(0,1);
        const std::complex<double> R = (Z_end - Z0)/(Z_end + Z0);
        for (size_t i=0; i<N; i++) {
            double xi = x[i];
                    
            // normalizing p(0) = 1
            // p[i] = (std::exp(-j*k*xi) + R*std::exp(-2.0*j*k*L)*std::exp(j*k*xi))/ (1.0 + R*std::exp(-2.0*j*k*L));
            // u[i] = (std::exp(-j*k*xi) - R*std::exp(-2.0*j*k*L)*std::exp(j*k*xi))/ (Z0*(1.0 + R*std::exp(-2.0*j*k*L)));
            
            //normalizing u(0) = 1

            p[i] = Z0 * (std::exp(-j*k*xi) + R*std::exp(-2.0*j*k*L)*std::exp(j*k*xi))/ (1.0 - R*std::exp(-2.0*j*k*L));
            u[i] = (std::exp(-j*k*xi) - R*std::exp(-2.0*j*k*L)*std::exp(j*k*xi))/ (1.0 - R*std::exp(-2.0*j*k*L));
        }

}

Geometry::TransferMatrix Cylinder::T(double x, double freq) const {
    Geometry::TransferMatrix M;
    double k = 2.0*constants::PI*freq/acoustics::c;
    const std::complex<double> j(0,1);
    M.t11 = std::cos(k*x);
    M.t12 = j*Z0*std::sin(k*x);
    M.t21 = j/Z0*std::sin(k*x);
    M.t22 = std::cos(k*x);
    return M;
}

void Cylinder::computeFieldFromMatrix(double freq, 
    std::complex<double> Z_end,
    std::vector<double>& x, 
    std::vector<std::complex<double>>& p, 
    std::vector<std::complex<double>>& u)
    const {
        size_t N = x.size();
        if (p.size() != N) {
            p.resize(N);
        }
        if (u.size() != N) {
            u.resize(N);
        }
        Geometry::TransferMatrix TL = T(L, freq);
        std::complex<double> denom = TL.t21*Z_end + TL.t22;
        for (size_t i=0; i<N; ++i) {
            double xi = x[i];
            // p[i] = (T(L-xi, freq).t11*Z_end + T(L-xi, freq).t12) / denom;
            // u[i] = (T(L-xi, freq).t21*Z_end + T(L-xi, freq).t22) / denom;

            Geometry::TransferMatrix T_xi = T(xi, freq);
            p[i] = (T_xi.t22*TL.t11*Z_end + T_xi.t22*TL.t12 - T_xi.t12*TL.t21*Z_end - T_xi.t12*TL.t22)/denom;
            u[i] = (T_xi.t11*TL.t21*Z_end + T_xi.t11*TL.t22 - T_xi.t21*TL.t11*Z_end - T_xi.t21*TL.t12)/denom;
        }
        
    }



//Cone
std::complex<double> Cone::inputImpedance(double freq, std::complex<double> Z_end) const {
    Geometry::TransferMatrix T_L = T(L, freq);
    return (T_L.t11*Z_end + T_L.t12) / (T_L.t21*Z_end + T_L.t22);
}

Geometry::TransferMatrix Cone::coneMatrix(double length, double r_start, double r_end, double freq) {
    Geometry::TransferMatrix M;
    double k = 2.0*constants::PI*freq/acoustics::c;
    const std::complex<double> j(0,1);

    if (length==0.0) {
        M.t11 = 1.0;
        M.t12 = 0.0;
        M.t21 = 0.0;
        M.t22 = 1.0;
    }
    else {
        double costh = std::sqrt(length*length - (r_end-r_start)*(r_end-r_start))/length;
        M.t11 = r_end/r_start*std::cos(k*length) - (r_end-r_start)/r_start*math::sinc(k*length);
        M.t12 = j*(1.0+costh)/2.0*acoustics::rho0*acoustics::c/(constants::PI*r_start*r_end)*std::sin(k*length);
        M.t21 = j*2.0/(1.0+costh)*constants::PI/(acoustics::rho0*acoustics::c) * (r_start*r_end*std::sin(k*length) + (r_end-r_start)*(r_end-r_start)*( math::sinc2(k*length) ) );
        M.t22 = r_start/r_end*std::cos(k*length) + (r_end-r_start)/r_end*math::sinc(k*length);
        // M.t21 = (M.t11*M.t22-1.0)/M.t12;
    }
    return M;
}

Geometry::TransferMatrix Cone::T(double x, double freq) const {
    double rx = R1 + (R2-R1)/L*x;
    return coneMatrix(x, R1, rx, freq);
}

void Cone::computeFieldFromMatrix(double freq, 
    std::complex<double> Z_end,
    std::vector<double>& x, 
    std::vector<std::complex<double>>& p, 
    std::vector<std::complex<double>>& u)
    const {
        size_t N = x.size();
        if (p.size() != N) {
            p.resize(N);
        }
        if (u.size() != N) {
            u.resize(N);
        }
        Geometry::TransferMatrix TL = T(L, freq); //T(L, R1, R2)
        std::complex<double> denom = TL.t21*Z_end + TL.t22;
        for (size_t i=0; i<N; ++i) {
            double xi = x[i];
            double rx = R1 + (R2-R1)/L*xi;

            // Geometry::TransferMatrix T_left = coneMatrix(L-xi, rx, R2, freq); //T(L-x[i]; rx, R2)
            // p[i] = (T_left.t11*Z_end + T_left.t12) / denom;
            // u[i] = (T_left.t21*Z_end + T_left.t22) / denom;


            // p[i] = (T(xi, freq).t22*T(L, freq).t11*Z_end + T(xi, freq).t22*T(L, freq).t12 - T(xi, freq).t12*T(L, freq).t21*Z_end - T(xi, freq).t12*T(L, freq).t22)/denom;
            // u[i] = (T(xi, freq).t11*T(L, freq).t21*Z_end + T(xi, freq).t11*T(L, freq).t22 - T(xi, freq).t21*T(L, freq).t11*Z_end - T(xi, freq).t21*T(L, freq).t12)/denom;


            Geometry::TransferMatrix T_xi = T(xi, freq);
            p[i] = (T_xi.t22*TL.t11*Z_end + T_xi.t22*TL.t12 - T_xi.t12*TL.t21*Z_end - T_xi.t12*TL.t22)/denom;
            u[i] = (T_xi.t11*TL.t21*Z_end + T_xi.t11*TL.t22 - T_xi.t21*TL.t11*Z_end - T_xi.t21*TL.t12)/denom;

        }
        
    }
