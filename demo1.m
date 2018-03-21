% demo1 calculate axial and lateral PSFs of axicon-based Bessel beam
% module with four differnt masks with the inner and the outer diameters 
% located at 1/e, 1/(5e), 1/(10e) and 1/(15e) of the peak amplitude of the
% ring, respectively.
%the results are saved in the subfolder (named "result") of the folder where codes are located.
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
offset_Lens2=0; % mm, offset of position of the lens 2; 0: the mask is at the focal plane of the Lens 2. If it is > 0, e.g., 20mm, Lens 2 is moved 20mm further away from the mask. 

% get the parameters for 4 masks
outputMask=maskDesign(beamD,f1,axicon,wavelength);
p=length(outputMask.outerD);
for jj=1:p
    if jj==1
        PSFz=cell(p,1);
        output=cell(p,1);
    end
    mask.outerD=outputMask.outerD(jj);
    mask.innerD=outputMask.innerD(jj);
   [PSFz{jj},output{jj}]=PSFofBesselBeam_Axicon(beamD,f1,mp,wavelength,x,y,z,axicon,obj,mask,offset_Lens2);
end
filePath=mfilename('fullpath');    
filePath=fileparts(filePath);
result='result';
resultPath=fullfile(filePath,result);
if exist(resultPath,'dir')
else
    mkdir(resultPath)
end        
save(fullfile(resultPath,'demo1_output.mat'));



    

