# Acoustical characterization of different geometries
Browser-based simulator for variable cylindrical and conical waveguides to explore how simple waveguide geometries influence acoustic resonances and impedance spectra. Backend written in C++ and compiled to WebAssembly via emscripten. (Runs fully client side.)

# Build
``` 
emcc geometries.cpp wasm_interface.cpp   
-s WASM=1  \
-s MODULARIZE=1   \
-s EXPORT_NAME='createModule' \ 
-s EXPORTED_FUNCTIONS="['_malloc','_free','_computeFieldCylinder', '_computeFieldCone', '_computeImpedanceAbsCylinder', '_computeImpedanceAbsCone', '_computeImpedanceSpectrumCylinder', '_computeImpedanceSpectrumCone']" \
-s EXPORTED_RUNTIME_METHODS="['cwrap','HEAPF64']" \
-s ALLOW_MEMORY_GROWTH=1 \
-flto \ 
-O3 \
-o ../build/acoustics.js
```

# Theory
Acoustic wave propagation is described by the Euler equation $\rho_0 \partial u/\partial t = -\partial p/\partial x$ with pressure $p$, particle velocity $u$ and equilibrium density $\rho_0$ as well as the linearized continuity equation $\partial \rho'/\partial t +\rho_0 \partial u / \partial x=0$ where $\rho'$ is the density perturbation. Assuming harmonic time dependence, together with the state equation $p=c^2 \rho'$, with the speed of sound $c$, this leads to the linear system  
$$
% \begin{align*}
% \frac{\partial \hat{p}}{\partial x} &= i \omega \frac{\rho_0}{S} \hat{U} \\
% \frac{\partial \hat{U}}{\partial x} &= -i \omega \frac{S}{\rho_0 c^2} \hat{p}
% \end{align*}

\begin{bmatrix}\partial \hat{p}/\partial x \\ \partial \hat{U}/\partial x\end{bmatrix} = 
\begin{bmatrix}0 & -i \omega \rho_0 /S \\-i \omega S/(\rho_0 c^2)  &0 \end{bmatrix} 
\begin{bmatrix}\hat{p} \\ \hat{U}\end{bmatrix}

$$
where $S$ is the (potentially variable) cross-sectional area of the waveguide and $U=Su$ is the volume velocity. $\hat{p}$ and $\hat{U}$ denote the respective complex amplitudes.

The linearity of the governing equations allows the acoustic state to be propagated using a transfer matrix T via

$$
\mathbf y(x_1) = T \mathbf y(x_2),
$$
where $\mathbf y(x) = [\hat{p}, \hat{U}]^\text{T}$ is the state vector. 
This is in principle possible for every geometry. For the simplest case where $S$ is constant, the transfer matrix can be calculated directly by evaluating the corresponding matrix exponential. For the case of a cylindrical waveguide, this yields 
$$
T_\text{cyl} = \begin{bmatrix} \cos(k L) & i \rho_0 c/(\pi R^2) \sin(kL) \\ i \pi R^2/(\rho_0 c) \sin(k L) & \cos(kL)\end{bmatrix}
$$
where R and L denote the cylinder radius and length, respectively, and $k=\omega/c$ is the wave number. For a cone, the transfer matrix is
$$
T_\text{cone} = 
\begin{bmatrix} 
\frac{R_2}{R_1}\cos(kL) - \frac{R_2-R_1}{R_1}\sin(kL)/(kL) & 
i\frac{(1+\cos\vartheta)}{2} \rho_0 c/(\pi R_1 R_2)\sin(kL) \\
i\frac{2}{1+\cos \vartheta}\pi/(\rho_0 c) \biggl(R_1 R_2 \sin(kL) + (R_2-R_1)^2\Big( \sin(kL)/(k^2L^2) - \cos(kL)/(kL)\Big) \biggl) & 
\frac{R_2}{R_1}\cos(kL) + \frac{R_2-R_1}{R_2}\sin(kL)/(kL)
\end{bmatrix} 
$$
(see A. Ernoult and J. Kergomard: [Acta Acustica **4**, 2 (2020)](https://doi.org/10.1051/aacus/2020005)). 


This formalism can be used to directly calculate the input impedance 
$$
Z_\mathrm{in} = \frac{\hat p_\mathrm{in}}{\hat U_\mathrm{in}},
$$
relating pressure and volume velocity at the input of a waveguide, as a measure for the general acoustic behavior of the geometry. Resonances of the acoustic system can be identified from extrema of $|Z_\mathrm{in}|$. Whether they are maxima or minima depends on the effective boundary condition at the observation point. If the input behaves approximately like a rigid termination (closed end), then $U_\mathrm{in} \approx 0$ at resonance. This corresponds to $|Z_\mathrm{in}|\to \infty$, so resonances appear as impedance peaks. If the input behaves approximately like an open end, then $p_\mathrm{in} \approx 0$ at resonance. This corresponds to $|Z_\mathrm{in}|\to 0$, so resonances appear as impedance minima.
