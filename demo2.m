% demo1 calculate axial and lateral PSFs of axicon-based Bessel beam
% module with the displacement of the lens L2 varying from -20 mm to 20 mm;
% the mask has the inner and the outer diameters 
% located at 1/e of the peak amplitude of the ring
%the results are saved in the subfolder (named "result") of the 
% folder where codes are located.
clear;
beamD=2.8;% mm, diameter of the beam at the axicon
f1=80;% mm, focal length of Lens L1 that is right after the axicon
mp=200/35.2*200/150;% magnification from the mask to the back focal plane of the objective
wavelength=0.94;%um,
x=0;
y=0;
z=-10:1:120; %um, range of z of the PSF relative to the focal plane
axicon.alpha=1; % alpha angle of axicon = (180-apex angle)/2
axicon.refind=1.4512;% refractive index of axicon material at the used wavelength; at 0.94
%     axicon.refind=1.4469;% refractive index of axicon material at the used wavelength;
axicon.diameter=25;% mm
% options: https://www.asphericon.com/wp-content/uploads/SPA-North-America_2017.pdf
%: 0.5, 1, 2, 5, 10, 20
obj.NA=0.95;% NA of objective
obj.magnification=25;% magnification magnification of objective
obj.refind=1.333;% refractive index of used medium for objective  
obj.tubeLength=200;% mm
offset_Lens2=-20:5:20; % mm, offset of position of the lens 2; 0: the mask is at the focal plane of the Lens 2. If it is > 0, e.g., 20mm, Lens 2 is moved 20mm further away from the mask. 

% calculate the inner and outer diameters that are located at 1/e of the peak amplitude of the ring
outputMask=maskDesign(beamD,f1,axicon,wavelength);
mask.outerD=outputMask.outerD(1);
mask.innerD=outputMask.innerD(1);
p=length(offset_Lens2);
PSFz=cell(p,1);
output=cell(p,1);    
for jj=1:p
   [PSFz{jj},output{jj}]=PSFofBesselBeam_Axicon(beamD,f1,mp,wavelength,x,y,z,axicon,obj,mask,offset_Lens2(jj));
end
filePath=mfilename('fullpath');    
filePath=fileparts(filePath);
result='result';
resultPath=fullfile(filePath,result);
if exist(resultPath,'dir')
else
    mkdir(resultPath)
end        
save(fullfile(resultPath,'demo2_output.mat'));



    

