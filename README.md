# simulation-for-axicon-based-Bessel-beam-module
This repository is to share the simulation codes for the paper “50 Hz volumetric functional imaging with continuously adjustable depth of focus” by Lu et al.  
Codes were written and tested in MATLAB 2016b. 
Supplementary Codes contain four codes: PSFofBesselBeam_Axicon.m, maskDesign.m, demo1.m, and demo2.m. 
The code “maskDesign.m” generates four masks with the inner and the outer diameters located at 1/e, 1/(5e), 1/(10e) and 1/(15e) of the peak amplitude of the ring, respectively. The transmittance ratios of these four masks under ideal conditions are 92.4%, 99.0%, 99.6% and 99.8%, respectively.
The code “PSFofBesselBeam_Axicon.m” simulate the electrical fields at the mask plane and at the back focal plane of the objective, axial PSF, PSF along y, and PSF along x.
The code “demo1.m” first generates 4 masks using the code “maskDesign.m” and then simulate corresponding results using the code “PSFofBesselBeam_Axicon.m”.
The code “demo2.m” calculate the PSFs with different displacements of lens L2 (refer to Figure 1 and Visualization 1 in the manuscript) using the code “PSFofBesselBeam_Axicon.m”.
