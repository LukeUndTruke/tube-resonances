# Acoustical characterization of different geometries
Browser-based simulator for variable cylindrical and conical waveguides to explore how simple waveguide geometries influence acoustic resonances and impedance spectra. Backend written in C++ and compiled to WebAssembly via emscripten. (Runs fully client side.) 

Live Demo: [https://lukeundtruke.github.io/tube-resonances](https://lukeundtruke.github.io/tube-resonances)


# Theory
Acoustic wave propagation is described by the Euler equation 

$$
\rho_0 \frac{\partial u}{\partial t} = -\frac{\partial p}{\partial x}
$$ 

with pressure $p$, particle velocity $u$ and equilibrium density $\rho_0$ as well as the linearized continuity equation 

$$\frac{\partial \rho'}{\partial t} +\rho_0 \frac{\partial u}{\partial x}=0
$$ 

where $\rho'$ is the density perturbation. Assuming harmonic time dependence, together with the state equation $p=c^2 \rho'$, with the speed of sound $c$, this leads to the linear system  

$$
\begin{bmatrix}\partial \hat{p}/\partial x \\\ \partial \hat{U}/\partial x\end{bmatrix} = 
\begin{bmatrix}0 & -i \omega \frac{\rho_0}{S} \\\ -i \omega \frac{S}{\rho_0 c^2}  &0 \end{bmatrix} 
\begin{bmatrix}\hat{p} \\\ \hat{U}\end{bmatrix}
$$

where $S$ is the (potentially variable) cross-sectional area of the waveguide and $U=Su$ is the volume velocity. $\hat{p}$ and $\hat{U}$ denote the respective complex amplitudes.

The linearity of the governing equations allows the acoustic state to be propagated using a transfer matrix T via

$$
\mathbf y(x_1) = T \mathbf y(x_2),
$$

where $\mathbf y(x) = [\hat{p}, \hat{U}]^\text{T}$ is the state vector. 
This is in principle possible for every geometry. For the simplest case where $S$ is constant, the transfer matrix can be calculated directly by evaluating the corresponding matrix exponential. For the case of a cylindrical waveguide, this yields 

$$
T_\text{cyl} = \begin{bmatrix} \cos(k L) & i \frac{\rho_0 c}{\pi R^2} {\sin(kL)} \\\ 
i \frac{\pi R^2}{\rho_0 c} \sin(k L) & \cos(kL)\end{bmatrix}
$$

where R and L denote the cylinder radius and length, respectively, and $k=\omega/c$ is the wave number. For a cone, the transfer matrix is

$$
T_\text{cone} = \begin{bmatrix} 
\frac{R_2}{R_1}\cos(kL) - \frac{R_2-R_1}{R_1}\frac{\sin(kL)}{kL} & 
i\frac{(1+\cos\vartheta)}{2} \frac{\rho_0 c}{\pi R_1 R_2}\sin(kL) \\\
i\frac{2}{1+\cos \vartheta}\frac{\pi}{\rho_0 c} \biggl(R_1 R_2 \sin(kL) + (R_2-R_1)^2\Big( \frac{\sin(kL)}{k^2L^2} - \frac{\cos(kL)}{kL}\Big) \biggl) & 
\frac{R_2}{R_1}\cos(kL) + \frac{R_2-R_1}{R_2}\frac{\sin(kL)}{kL}
\end{bmatrix} 
$$

where $\vartheta$ is half the cone's opening angle (see A. Ernoult and J. Kergomard: [Acta Acustica **4**, 2 (2020)](https://doi.org/10.1051/aacus/2020005)). 


This formalism can be used to directly calculate the input impedance 

$$
Z_\mathrm{in} = \frac{\hat p_\mathrm{in}}{\hat U_\mathrm{in}},
$$

relating pressure and volume velocity at the input of a waveguide, as a measure for the general acoustic behavior of the geometry. Resonances of the acoustic system can be identified from extrema of $|Z_\mathrm{in}|$. Whether they are maxima or minima depends on the effective boundary condition at the observation point. If the input behaves approximately like a rigid termination (closed end), then $U_\mathrm{in} \approx 0$ at resonance. This corresponds to $|Z_\mathrm{in}|\to \infty$, so resonances appear as impedance peaks. If the input behaves approximately like an open end, then $p_\mathrm{in} \approx 0$ at resonance. This corresponds to $|Z_\mathrm{in}|\to 0$, so resonances appear as impedance minima.

# Interesting physical phenomena
### Influence of conicity
As stated above, the natural resonant frequencies of the geometry with a closed input (left) end are the input impedance peaks, and the natural frequencies of the geometry with an open end are the input impedance minima. For more detail regarding this following discussion, see R. D. Ayers, L. J. Eliason, and D. Mahgerefteh: [American Journal of Physics **53**, 6: 528–537 (1985)](https://doi.org/10.1119/1.14233). 

#### Open right end ($Z_r=0$)
Introducing a conicity to the geometry by changing $R_1$ or $R_2$ such that $R_1<R_2$ shows that, for zero termination impedance $Z_r=0$ (open right end), the resonances of the structure with an open input (left) end do not depend on the conicity; the resonances of an open-open cone are the same as those of an open-open cylinder of the same length. 

The resonances of the structure with a closed left end, however, do depend on conicity. For $R_1=R_2$, they are of course the same as those of a cylinder closed at its left end. With decreasing $R_1$, they tend towards those of the structure with an open left end. The limiting case, $R_1 \to 0$, corresponds to a complete cone. In this case, the distinction between an open and a closed left end becomes meaningless, and the resonant frequencies of the structure with an open and with a closed left end coincide.

#### Closed right end ($Z_r \to \infty$)
In this case, introducing a finite conicity to the cylinder ($R_1<R_2$) does not affect the higher modes of the structure with an open left end, However, it does affect the very first mode. Its frequency drops to $0$ when the cone is closed completely ($R_1 \to 0$), such that it disappears completely. This can be understood by looking at the shape of the standing wave of this resonance. For this resonance, there is only one velocity antinode (pressure node) in the whole geometry, at the left end. When closing this end by reducing $R_1 \to 0$, this velocity antinode cannot exist anymore, so the resonance disappears. The higher modes are not affected, because for these resonances, there is at least one node/antinode _inside_ the geometry. The boundary condition on the left is of course distorted by introducing a conicity; for a complete cone, the pressure node/velocity antinode at the left end turns into a pressure antinode/velocity node. However, these resonances can exist. 

For a structure that is closed also at its left end (i.e., a completely closed cylinder/cone), the conicity does influence the resonance frequencies. Just as in the case $Z_r=0$, they are shifted towards those of the structure with an open left end, and in the limiting case $R \to 0$ (complete cone), the resonances of the structure with an open and with a closed left end coincide.

### Influence of termination (load) impedance
Zero or infinite real part $Z_r$ of termination impedance means an open or closed end, respectively. For a cylinder, the reflection coefficient $R$ (the ratio between reflected and incident pressure amplitude) is 

$$
\Gamma = \frac{Z_L - Z_0}{Z_L + Z_0}.
$$

$Z_0$ is the so-called characterisic impedance of the medium inside the geometry, defined by $Z_0 = \hat{p}/\hat{U}$ for a pure travelling wave, yielding $Z_0 = \rho_0 c / (\pi R^2)$ for a cylinder. (This can be derived from the above equations.) This means that for both an open and a closed right end, $|\Gamma|=1$, the incident wave is fully reflected, leading to a pure standing wave. Chosing $Z_L = Z_0$ means that the load/termination impedance of the material to the right of the geometry is the same as that inside the geometry. The result is then a pure travelling wave, since the incident wave simply continues to propagate to the right into the same medium. 


# Build

``` 
emcc geometries.cpp wasm_interface.cpp \   
-s WASM=1 \
-s MODULARIZE=1 \
-s EXPORT_NAME='createModule' \ 
-s EXPORTED_FUNCTIONS="['_malloc','_free','_computeFieldCylinder', '_computeFieldCone', '_computeImpedanceAbsCylinder', '_computeImpedanceAbsCone', '_computeImpedanceSpectrumCylinder', '_computeImpedanceSpectrumCone']" \
-s EXPORTED_RUNTIME_METHODS="['cwrap','HEAPF64']" \
-s ALLOW_MEMORY_GROWTH=1 \
-flto \ 
-O3 \
-o ../build/acoustics.js
```