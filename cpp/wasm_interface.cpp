#include <vector>
#include <complex>
#include "geometries.h"

extern "C"{
    void computeFieldCylinder(double L, double R, double freq, double Zr, double Zi, 
        double* x, int N, double* p_abs, double* p_arg, double* u_abs, double* u_arg){
            Cylinder cyl(L, R);
            
            std::complex<double> Z_end(Zr, Zi);
            std::vector<double> x_vec(x, x+N);
            std::vector<std::complex<double>> p_vec(N);
            std::vector<std::complex<double>> u_vec(N);

            cyl.computeFieldFromMatrix(freq, Z_end, x_vec, p_vec, u_vec);

            for (int i=0; i<N; ++i){
                p_abs[i] = std::abs(p_vec[i]);
                p_arg[i] = std::arg(p_vec[i]);
                u_abs[i] = std::abs(u_vec[i]);
                u_arg[i] = std::arg(u_vec[i]);
            }
         }
    
    void computeFieldCone(double L, double R1, double R2, double freq, double Zr, double Zi,
        double* x, int N, double* p_abs, double* p_arg, double*u_abs, double* u_arg){
            Cone cone(L, R1, R2);
            
            std::complex<double> Z_end(Zr, Zi);
            std::vector<double> x_vec(x, x+N);
            std::vector<std::complex<double>> p_vec(N);
            std::vector<std::complex<double>> u_vec(N);

            cone.computeFieldFromMatrix(freq, Z_end, x_vec, p_vec, u_vec);

            for (int i=0; i<N; ++i){
                p_abs[i] = std::abs(p_vec[i]);
                p_arg[i] = std::arg(p_vec[i]);
                u_abs[i] = std::abs(u_vec[i]);
                u_arg[i] = std::arg(u_vec[i]);
            }
        }
    
    double computeImpedanceAbsCylinder(double L, double R, double freq, double Zr, double Zi) {
        Cylinder cyl(L, R);
        std::complex<double> Z_end(Zr, Zi);
        return std::abs(cyl.inputImpedance(freq, Z_end));
    }
    double computeImpedanceAbsCone(double L, double R1, double R2, double freq, double Zr, double Zi) {
        Cone cone(L, R1, R2);
        std::complex<double> Z_end(Zr, Zi);
        return std::abs(cone.inputImpedance(freq, Z_end));
    }

    void computeImpedanceSpectrumCylinder(double L, double R, double* freqArray, int NFreq, double Zr, double Zi,
        double* out_abs, double* out_arg){
        Cylinder cyl(L, R);

        std::complex<double> Z_end(Zr, Zi);
        std::complex<double> Z_in;

        for (size_t i=0; i<NFreq; i++) {
            Z_in =  cyl.inputImpedance(freqArray[i], Z_end);
            out_abs[i] = std::abs(Z_in);
            out_arg[i] = std::arg(Z_in);
        }
    }

    void computeImpedanceSpectrumCone(double L, double R1, double R2, double* freqArray, int NFreq, double Zr, double Zi,
        double* out_abs, double* out_arg){
        Cone cone(L, R1, R2);

        std::complex<double> Z_end(Zr, Zi);
        std::complex<double> Z_in;

        for (size_t i=0; i<NFreq; i++) {
            Z_in =  cone.inputImpedance(freqArray[i], Z_end);
            out_abs[i] = std::abs(Z_in);
            out_arg[i] = std::arg(Z_in);
        }
    };
}