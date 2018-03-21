function [PSFz,output]=PSFofBesselBeam_Axicon(beamD,f1,mp,wavelength,x,y,z,axicon,obj,mask,offset_Lens2)
% PSFofBesselBeam_Axicon calculates PSF near focal point using Richards and
% Wolf integrals. 
% the optical setup: axicon-->Lens 1 -->mask-->relay lenses L2 and L3
% -->scanners -->objective -->samples;
% input argument beamD is the diameter of the Gaussian beam at axicon;
% f1 is the focal length of lens after the axicon;
% mp is the magnification from the mask plane (Fourier plane of axicon) to back
% pupil plane of objective;
% wavelength is the two-photon excitation wavelength;
% x,y,z are coordinates relative to focal point in um
% axicon includes the paraemters of axicon
% axicon.alpha is alpha angle of axicon which is 0.5*(180-apexAngle)
% axicon.refind is the refractive index of axicon material
% axicon.diameter is the diameter of the axicon
% obj includes paramters of obj
%   obj.tubeLength tubelength of objective in mm - Nikon objectives: 200 mm;
%   Olympus objectives: 180 mm; Zeiss objective: 165 mm. 
%   obj.NA NA of objective
%   obj.magnification magnification of objective
%   obj.refind, refractive index of immersion medium used for objective
%
% mask includes parameters of annular filter
%	mask.outerD outer diameter of annular ring in mm
%   mask.innerD inner diameter of annular ring in mm
%   The output arguments contain PSF and output
%
%-------------example------------------
%     beamD=2.8;% mm, diameter of the beam at the axicon
%     f1=80;% mm, focal length of Lens L1 that is right after the axicon
%     mp=200/35.2*200/150;% magnification from the mask to the back focal plane of the objective
%     wavelength=0.94;%um,
% 
%     x=0;
%     y=0;
% %     z=-120:0.5:10;
%     z=-10:1:150; %um, range of z of the PSF relative to the focal plane
%     axicon.alpha=1; % alpha angle of axicon;
%     axicon.refind=1.4512;% refractive index of axicon material at the used wavelength; at 0.94
% %     axicon.refind=1.4508; % at 970nm
%     axicon.diameter=25;% mm
% % options: https://www.asphericon.com/wp-content/uploads/SPA-North-America_2017.pdf
% %: 0.5, 1, 2, 5, 10, 20
%     obj.NA=0.95;% NA of objective
%     obj.magnification=25;% magnification magnification of objective
%     obj.refind=1.333;% refractive index of used medium for objective  
%     obj.tubeLength=200;% mm
%     mask.outerD=1.317;% mm, outer diameter of annular ring
%     mask.innerD=1.203;% mm, inner diameter of annular ring    
%     offset_Lens2=-12; % mm, displacement of lens L2
%     [PSFz,output]=PSFofBesselBeam_Axicon(beamD,f1,mp,wavelength,x,y,z,axicon,obj,mask,offset_Lens2);
    
if nargin==0
    beamD=2.8;% mm, diameter of the beam at the axicon
    f1=80;% mm, focal length of Lens L1 that is right after the axicon
    mp=200/35.2*200/150;% magnification from the mask to the back focal plane of the objective
    wavelength=0.94;%um,

    x=0;
    y=0;
%     z=-120:0.5:10;
    z=-10:1:150; %um, range of z of the PSF relative to the focal plane
    axicon.alpha=1; % alpha angle of axicon;
    axicon.refind=1.4512;% refractive index of axicon material at the used wavelength; at 0.94
%     axicon.refind=1.4508; % at 970nm
    axicon.diameter=25;% mm
% options: https://www.asphericon.com/wp-content/uploads/SPA-North-America_2017.pdf
%: 0.5, 1, 2, 5, 10, 20
    obj.NA=0.95;% NA of objective
    obj.magnification=25;% magnification magnification of objective
    obj.refind=1.333;% refractive index of used medium for objective  
    obj.tubeLength=200;% mm
    mask.outerD=1.317;% mm, outer diameter of annular ring
    mask.innerD=1.203;% mm, inner diameter of annular ring    
    offset_Lens2=-12; % mm, offset of position of the lens 2; 0: the mask is at the focal plane of the Lens 1. If it is > 0, e.g., 20mm, Lens 2 is moved 20mm further away from the mask. 
end
%% save input information
output.mask=mask;
output.f1=f1;
output.beamD=beamD;
output.wavelength=wavelength;
output.mp=mp;
output.axicon=axicon;
output.obj=obj;
output.offset_Lens2=offset_Lens2;
output.x=x;
output.y=y;
output.z=z;
%% simulation parameters at the plane of the axicon
step_axicon=0.4;% um, step size of each pixel for the coordinate on the axicon
pixelNumber=axicon.diameter*1000/step_axicon;% pixel number for the coordinate on the axicon
halfPixelNumber=pixelNumber*1/2;% pxiel numbers from center of the axicon to the edge
beta=(axicon.refind-1)*axicon.alpha/180*pi;% deflected angle of the beam

%% objective parameters    
% tube lens:Nikon 200 mm, Olympus: 180mm, Zeiss: 165 mm;
tubeLength=obj.tubeLength;
NA=obj.NA;
magnification=obj.magnification;
refind=obj.refind;
fOB=tubeLength/magnification; % focal length of the objective
rOB=tubeLength/magnification*NA; % radius of objective back pupil
%% mask parameters
outerD=mask.outerD; 
innerD=mask.innerD;
r2=outerD/2;
r1=innerD/2;
% coordinate at mask
step_mask=0.4;% micron;
r_mask=0:step_mask:rOB/mp*1000*2; % radial coordinate at mask
initialPhase=0; % 
beamr=beamD/2*1000; % obtain beam radius at the axicon in the unit of um
r_axicon=(0:1:halfPixelNumber-1)*step_axicon;% radial coordinate at axicon

%% field at the axicon
% note: for the convex axicon there is a minus sign in the equation. we can
% refer to the function of the lens analogously which is
% exp(-1i*2*pi*r^2/2/f), where f is the focal length of the lens. 
field_axicon=exp(-1i*(2*pi*r_axicon*beta/wavelength-initialPhase));
field_axicon=field_axicon.*exp(-r_axicon.^2/beamr/beamr);% phase added upon the Gaussian beam
%% field at the mask
field_mask=fourier_circular(field_axicon,wavelength,r_axicon,f1,r_mask);
field_mask_save=field_mask;
%% calculate the effective NA
ampRatio=1/exp(1);
ampRatioX=findXatAmpRatio(abs(field_mask),ampRatio);
NA_effective=r_mask(ampRatioX)/1000*mp/fOB;
if NA_effective>obj.NA
    warndlg(['NA of Bessel is',num2str(NA_effective,'%.2f'), ' bigger than NA of obj (',num2str(obj.NA,'%.2f'),')'],' Warning ')
end
output.NA_effective=NA_effective;
%% calculate Transmittance after the mask
powerTmp=abs(field_mask).^2.*r_mask;
flag1=r_mask*mp<=rOB*1000;%coordinates that are inside the objective  
flag2=(r_mask>=innerD/2*1000)&(r_mask<=outerD/2*1000); % coordinates that are inside the annular filter    
residualPower_mask=sum(powerTmp(flag2))/sum(powerTmp);
%% plot electrcial field at mask (Fourier plane of axicon)
% limab=[0,f1*beta*2.5];
limab=[0,rOB/mp*1.2];
IntensityY=abs(field_mask_save);maxY=max(IntensityY);
IntensityY=IntensityY/maxY;
% IntensityY=IntensityY.^2;
figure(11);clf;
panelM=4;panelN=2;
set(gcf,'position',[27    77   839   592])


lineW=2;
titleSize=10;
labelSize=10;
plotOrder=1:panelM*panelN;
plotOrder=reshape(plotOrder,panelN,panelM);
plotOrder=plotOrder.';
plotOrder=plotOrder(:);
plot_i=1;
subplot(panelM,panelN,plotOrder(plot_i));
plot_i=plot_i+1;
plot(r_mask/1000,IntensityY,'b','LineWidth',lineW);
xlim(limab);
xlabel(['Transmittance=',num2str(residualPower_mask*100,'%.1f'),'%'])

h_title=title('Amp. & phase of E at Fourier plane of axicon');h_title.FontSize=titleSize;
% set(gca,'xTick',[0:0.3:1.8])
set(gca,'FontName','Arial')
b=1.2;
ylim([0,1.2]);set(gca,'yTick',[0,0.5,1])
set(gca,'color','none')
set(gca,'TickDir','out');
box off;
hold on; 
plot([r1,r1],[0,1]*b,'k--','LineWidth',lineW-0.5);
plot([r2,r2],[0,1]*b,'k--','LineWidth',lineW-0.5)
angleY=angle(field_mask_save)/pi*180;
angleY(abs(angleY-360+180)<=10^-6)=-180;

ylabel('Amplitude','Fontname','Arial','Fontsize',labelSize)
% xlabel('mm','Fontname','Arial','Fontsize',11);
xlim(limab);
subplot(panelM,panelN,plotOrder(plot_i));
plot_i=plot_i+1;
plot(r_mask/1000,angleY,'b','LineWidth',lineW);
% set(gca,'xTick',[0:0.3:1.8])
ylim([-200,200])
b=200;
hold on; 
plot([r1,r1],[-1,1]*b,'k--','LineWidth',lineW-0.5);
hLine1=plot([r2,r2],[-1,1]*b,'k--','LineWidth',lineW-0.5);
hLegend=legend([hLine1],'Inner and outer edges of annular mask');
hLegend.Position=[0.25 0.6956 0.1756 0.0212];
box off;
set(gca,'FontName','Arial')
set(gca,'color','none')
set(gca,'TickDir','out');
set(gcf,'PaperPositionMode','auto')
set(gca,'yTick',[-180,0,180])
ylabel('Phase(Degree)','Fontname','Arial','Fontsize',labelSize)
xlim(limab);
xlabel('r(mm)','Fontname','Arial','Fontsize',labelSize);


field_mask(~flag2)=0; % electrical field right after the mask;
%% adding defocus at the backpupil plane
r_mask=r_mask*mp; % change radial coordinate from plane of mask to objective pupil plane
offset_objective=offset_Lens2*mp^2; % unit mm; convert offset of lens 2 to 
% the offset of ring with respect to the objective; 
% Let us assume the Lens 2 is moved 20 mm away from the mask. It means
% distance between the mask and the Lens 2 is increased by 20 mm.
% Therefore, the conjugated plane after the tube lens should be 
% 20 mm*mp^2 closer to the Tube lens. To get the electrical field at the
% back pupil plane, we need to let electrical field proporgate with the
% distance of 20 mm*mp^2. 
field_pupil=proporgation_circularSymmetry(field_mask,wavelength,r_mask,offset_objective*1000);
%Note: to speed up the simulation, the distance between the Lens 2 and Lens
%3 is assumed to be fixed, which only displaces the PSF in few
%microns along axial direction. 

%% Transmittance after cut by the pupil edge
powerTmp=abs(field_pupil).^2.*r_mask;
residualPower_pupil=sum(powerTmp(flag1))/sum(powerTmp);

field_pupil_amp=abs(field_pupil(:));
field_pupil_phase=angle(field_pupil(:))/pi*180;

%% plot amplitude & phase at the pupil plane of the objective
subplot(panelM,panelN,plotOrder(plot_i));
plot_i=plot_i+1;
plot(r_mask/1000,field_pupil_amp/max(field_pupil_amp),'b','LineWidth',lineW);

xlabel(['Transmittance=',num2str(residualPower_pupil*100,'%.1f'),'%'])
h_title=title(['Amp. & Phase of E at pupil plane, D=',num2str(offset_Lens2,'%.1f'),'mm']);h_title.FontSize=titleSize;

set(gca,'FontName','Arial')
b=1.2;
ylim([0,1.2]);set(gca,'yTick',[0,0.5,1])
set(gca,'color','none')
set(gca,'TickDir','out');
box off;
hold on; 
plot([1,1]*rOB,[0,1],'r-','LineWidth',lineW);
ylabel('Amplitude','Fontname','Arial','Fontsize',labelSize)
xlim(limab*mp);
subplot(panelM,panelN,plotOrder(plot_i));
plot_i=plot_i+1;
plot(r_mask/1000,field_pupil_phase,'b','LineWidth',lineW);
ylim([-200,260])
hold on; 
hLine2=plot([1,1]*rOB,[-200,200],'r-','LineWidth',lineW);
hLegend=legend(hLine2,'Edge of pupil');
hLegend.Position=[0.3836    0.2535    0.0793    0.0212];
box off;
set(gca,'FontName','Arial')
set(gca,'color','none')
set(gca,'TickDir','out');
set(gcf,'PaperPositionMode','auto')
set(gca,'yTick',[-180,0,180])
ylabel('Phase(Degree)','Fontname','Arial','Fontsize',labelSize)
xlim(limab*mp);
xlabel('r(mm)','Fontname','Arial','Fontsize',labelSize);

%% get PSF near the focal plane
r_angle=asin(r_mask/fOB/refind/1000);
r_angle(~flag1)=[];
field_pupil(~flag1)=[];
[PSF] = Calc_Annular_Field_Integrals(x, y, z,field_pupil,wavelength,r_angle,refind);

PSF=squeeze(PSF);
xlen=length(x);
ylen=length(y);
zlen=length(z);
% figure(2);clf;
xyzlen=sort([xlen,ylen,zlen]);

output.maxPSF=max(PSF);
PSFz=PSF;
if xlen==1 && ylen==1 && zlen>1
    %% plot axial PSF
    PSF2=PSF/max(PSF);
    drawnow
    plotTmp=subplot(panelM,panelN,[plotOrder(plot_i),plotOrder(plot_i+1)]);
    set(plotTmp,'Position',[0.5703    0.58    0.3347    0.334])
    plot_i=plot_i+2;
     plot(z,PSF2,'b','LineWidth',lineW);
    box off;
    set(gca,'FontName','Arial')
    set(gca,'color','none')
    set(gca,'TickDir','out');
    set(gcf,'PaperPositionMode','auto')
    set(gca,'yTick',[0:0.5:1])     
    ylim([0, 1.2])
     xlabel('z(\mum)','Fontname','Arial','Fontsize',labelSize);
     ylabel('Axial PSF','Fontname','Arial','Fontsize',labelSize)
     xlim([z(1),z(end)]);
     flag=PSF2>=0.5;
     II=find(flag);     
     FWHMz=z(II(end))-z(II(1));   
     h_title=title(['Axial PSF; Lens 2 offset=',num2str(offset_Lens2,'%.1f'), 'mm;FWHM=',num2str(FWHMz),'\mum']);h_title.FontSize=titleSize;
      drawnow; 
    output.PSFzRaw=PSF;
    output.FWHMz=FWHMz;
    %% calculating lateral FWHM along y
    [~,I]=max(PSF2);
    output.PSFz=[z(:),PSF2(:)];
    z=z(I(1));
    x=0;
    y=-2:.02:2;
    [PSFy] = Calc_Annular_Field_Integrals(x, y, z,field_pupil,wavelength,r_angle,refind);
    PSFy=squeeze(PSFy);
    PSFy2=PSFy/max(PSFy);
    % fitting
    y2=y;y2(y2==0)=eps;
    w=2*pi*obj.NA/wavelength;
    %% do fitting to get more accurate lateral FWHM
    fun=@(par,z,r)(z-(par(1)*besselj(1,par(2)*r)./r/2/pi).^4);
    fun2=@(par,r)((par(1)*besselj(1,par(2)*r)./r/2/pi).^4);
    par0=[1,1]*w;
    par=lsqnonlin(fun,par0,[],[],[],PSFy2(:),y2(:));   
    output.PSFy=[y(:),PSFy2(:)];
    step=0.001;
    y_finner=0:step:2;y_finner(1)=eps;
    y_finner=[-y_finner(end:-1:1),y_finner(2:end)];
    PSFy_fit=fun2(par,y_finner);
    flag=PSFy_fit>=0.5;
    FWHM2=sum(flag)*(y_finner(2)-y_finner(1))*1;  
    output.PSFy_fit=[y_finner(:),y_finner(:)];
    output.FWHMy_fit=FWHM2;
    drawnow
    subplot(panelM,panelN,plotOrder(plot_i));
    plot_i=plot_i+1;
    plot(y,PSFy2,'b-','LineWidth',lineW);
    h_title=title(['PSF along y, FWHM=',num2str(FWHM2,'%.3f'),'\mum']);h_title.FontSize=titleSize;
    box off;
    set(gca,'FontName','Arial')
    set(gca,'color','none')
    set(gca,'TickDir','out');
    set(gca,'yTick',[0:0.5:1])     
    ylim([0, 1.2])
     xlabel('y(\mum)','Fontname','Arial','Fontsize',labelSize);
     ylabel('Vertical PSF','Fontname','Arial','Fontsize',labelSize)
     xlim([y(1),y(end)]);
 %% calculating along x axis    
    y=0;
    x=-2:.02:2;
    [PSFx] = Calc_Annular_Field_Integrals(x, y, z,field_pupil,wavelength,r_angle,refind);
    PSFx=squeeze(PSFx);
    PSFx2=PSFx/max(PSFx);
    y2=x(:);
    y2(y2==0)=eps;
    par2=lsqnonlin(fun,par,[],[],[],PSFx2(:),y2(:));
    output.PSFx=[x(:),PSFx2(:)];
    x_finner=y_finner;
    PSFx_fit=fun2(par2,x_finner);
     % optional: convolution with 0.2 um bead
%      PSFx_fit=convolvedWithp2umBead(PSFx_fit,step);
    flag=PSFx_fit>=0.5;
    FWHM2=sum(flag)*(x_finner(2)-x_finner(1))*1;  
    output.PSFx_fit=[x_finner(:),x_finner(:)];
    output.FWHMx_fit=FWHM2;
    drawnow
    subplot(panelM,panelN,plotOrder(plot_i));
%     plot_i=plot_i+1;
    plot(x,PSFx2,'b-','LineWidth',lineW);
    h_title=title(['PSF along x, FWHM=',num2str(FWHM2,'%.3f'),'\mum']);h_title.FontSize=titleSize;
    box off;
    set(gca,'FontName','Arial')
    set(gca,'color','none')
    set(gca,'TickDir','out');
    set(gca,'yTick',[0:0.5:1])     
    ylim([0, 1.2])
     xlabel('x(\mum)','Fontname','Arial','Fontsize',labelSize);
     ylabel('Horizontal PSF','Fontname','Arial','Fontsize',labelSize)
     xlim([x(1),x(end)]);
     drawnow
    

    filePath=mfilename('fullpath');    
    filePath=fileparts(filePath);
    result='result';
    resultPath=fullfile(filePath,result);
    if exist(resultPath,'dir')
    else
        mkdir(resultPath)
    end
%% save figure
    file=['outerD',num2str(outerD,'%.3f'),'mm_innerD',num2str(innerD,'%.3f'),'mm_'];
    file=['f1_',num2str(f1),'mm_',file];
    if offset_Lens2<0
        file=[file,'D_minus',num2str(-offset_Lens2,'%04.1f'),'mm.png'];
    else
        file=[file,'D_',num2str(offset_Lens2,'%04.1f'),'mm.png'];
    end

    file2=fullfile(resultPath,file);
    figure(11);
    savefig(gcf,strrep(file2,'.png','.fig'));
    print('-r100',file2,'-dpng')

elseif xyzlen(1)==1 && xyzlen(2)>1
    plotTmp=subplot(panelM,panelN,[plotOrder(plot_i),plotOrder(plot_i+1)]);
    set(plotTmp,'Position',[0.5703    0.58    0.3347    0.334])
    plot_i=plot_i+2;    
    imshow(PSF,[])
    h_title=title('PSF');h_title.FontSize=titleSize;
    [m,n]=size(PSF);
elseif zlen==1
    if xlen==1
        t=y;
    elseif ylen==1
        t=x;
    end
    plotTmp=subplot(panelM,panelN,[plotOrder(plot_i),plotOrder(plot_i+1)]);
    set(plotTmp,'Position',[0.5703    0.58    0.3347    0.334])
    plot_i=plot_i+2;    
    plot(t,PSF)
    xlabel('z (\mum)');
    ylabel('two photon PSF, lateral plot')
    PSF2=PSF/max(PSF);
    flag=PSF2>=0.5;
    FWHM=sum(flag)*(t(2)-t(1));
    h_title=title(['Lateral PSF, FWHM=',num2str(FWHM),'\mum']);h_title.FontSize=titleSize;
end

function field_mask=fourier_circular(field_axicon,wavelength,r_axicon,f,r_mask)
f=f*1000;
p=length(r_mask);  
field_mask=zeros(size(r_mask));
for ii=1:p
    a=field_axicon.*besselj(0,2*pi*r_axicon*r_mask(ii)/wavelength/f)*2*pi.*r_axicon;
    field_mask(ii)=trapz(r_axicon,a);
end

function [PSF,I0, I1, I2] = Calc_Annular_Field_Integrals(x, y, z,mask_obj,wavelength,r_angle,refind)
%CALC_ANNULAR_FIELD_INTEGRALS  Calculates the 3 Richards & Wolf integrals used to calc. EM field
%   CALC_ANNULAR_FIELD_INTEGRALS calculates the three integrals I0, I1,
%      and I2 used by Richards and Wolf in the determination of all six
%      components of the EM field in a vicinity of a focus, for a single
%      condition of annular illumination
%input arguments x, y, and z are coordinates near focal point; mask_obj is
%the eletrical field at the back focal plane of objective, innerD and
%outerD define inner diameter and outer diameter of annular ring. 
% 
xlen = length(x);
ylen = length(y);
zlen = length(z);
%
%create cylindrical coord. dimensionless arrays:
xmat = (x.') * ones(1, ylen);   %xlen by ylen matrix, each column filled with x vector
ymat = ones(xlen, 1) * y;       %xlen by xlen matrix, each row filled with y vector
rmat = sqrt((xmat .^ 2) + (ymat .^ 2)) + eps;  %projection of vector from focus to field pt onto xy plane
k=2*pi/wavelength;
uvec = k .* i .* refind .* z;   %uvec = ikz, k = mag. of wavevector in air
vmat = k .* refind .* rmat;  %vmat = k * rperp


%%
%
%define matrices for the definite integrals I0, I1, and I2:
I0 = ones(xlen, ylen, zlen);
I1 = ones(xlen, ylen, zlen);
I2 = ones(xlen, ylen, zlen);
phi=ones(xlen, ylen, zlen);
%
%define functions for integrands of I0, I1, and I2:
I0integrand = inline('sqrt(cos(th)) .* sin(th) .* (1 + cos(th)) .* besselj(0, (vm .* sin(th))) .* exp(uv .* cos(th))', 'th', 'uv', 'vm');
I1integrand = inline('sqrt(cos(th)) .* (sin(th) .^ 2) .* besselj(1, (vm .* sin(th))) .* exp(uv .* cos(th))', 'th', 'uv', 'vm');
I2integrand = inline('sqrt(cos(th)) .* sin(th) .* (1 - cos(th)) .* besselj(2, (vm .* sin(th))) .* exp(uv .* cos(th))', 'th', 'uv', 'vm');
%
interr = 1e-3;  %fractional error permitted for integral evaluation using quad
nullv = [];     %null vector needed for trace argument in quad fcn
%
z_save=z;
for z = 1:zlen  %loop thru all z positions
%     disp(['z=',num2str(z_save(z))])
%    z
    for x = 1:xlen  %loop thru all x positions
%        x
        for y = 1:ylen  %loop thru all y positions
            a=I0integrand(r_angle,uvec(z),vmat(x, y)).*mask_obj;
            b=I1integrand(r_angle,uvec(z),vmat(x, y)).*mask_obj;
            c=I2integrand(r_angle,uvec(z),vmat(x, y)).*mask_obj;           
            I0_tmp=trapz(a);
            I1_tmp=trapz(b);
            I2_tmp=trapz(c);
            I0(x,y,z)=I0_tmp;
            I1(x,y,z)=I1_tmp;
            I2(x,y,z)=I2_tmp;    
            phi(x,y,z)=angle(x+1i*y);
        end
    end
end
% phi=0;
Ex=(I0+I2.*cos(2*phi));
Ey=(I2.*cos(2*phi));
Ez=(-2i*I1.*cos(phi));
PSF=abs(Ex).^2+abs(Ey).^2+abs(Ez).^2;
PSF=PSF.^2;
function PSF=convolvedWithp2umBead(PSF,step)
beadSize=0.2;
siz=beadSize/step;
siz_y=1:siz;siz_y=siz_y-0.5*siz;
siz_y=siz_y*step;
r=beadSize/2;
kernel=(r^2-siz_y.^2).^0.5;
kernel=kernel/sum(kernel);
PSF=imfilter(PSF(:),kernel(:),'replicate');
PSF=PSF/max(PSF);
function [ampRatioX,halfThickness]=findXatAmpRatio(amplitude,ampRatio)

    m=length(ampRatio);
    ampRatioX=zeros(size(ampRatio));
    [peakV,peakI]=max(amplitude);
    if length(peakI)>1
        peakI=peakI(1);
        peakV=peakV(1);
    end
    amplitudeTmp=amplitude(peakI:end);
    for ii=1:m
        ampRatioValue=peakV*ampRatio(ii);
        amplitudeTmp2=abs(amplitudeTmp-ampRatioValue);
        [~,tmp]=min(amplitudeTmp2);
        ampRatioX(ii)=tmp(1);
        
    end
    ampRatioX=ampRatioX+peakI-1;
    halfThickness=ampRatioX-peakI;
    
function maskZ=proporgation_circularSymmetry(mask,wavelength,r,z)
% PROPORGATION_CIRCULARSYMEMETRY is to proporgate the electrical field mask
% by z distance. wavelength, r and z are in units of microns. The
% electrical field is circularly symmetric. r is the vector with the same
% size of mask, and it is the radius in polar coordinates. 
% alpha=linspace(0,15,2*10^4);% 0-15 degrees, sampled at 15/2/10^4 degree;
alpha=linspace(0,5,10^4);% 0-5 degrees, sampled at 5/10^4 degree;
alpha=alpha/180*pi;% convert degrees to radians

rho=sin(alpha)/wavelength;% cycles/micron, coordinate in Fourier plane;
p=length(rho);  
if z==0
    % to speed up
    maskZ=mask;
else
    %% Fourier-Bessel transform, get Fz0, angular spectra at z0
    Fz0=zeros(size(rho));
    for ii=1:p
        a=mask.*besselj(0,2*pi*r*rho(ii))*2*pi.*r;
        Fz0(ii)=trapz(r,a);
    end
    %% proporgation, get Fz, angular spectra at z
    Fz=Fz0.*exp(1i*2*pi*z.*sqrt(1-wavelength^2*rho.^2));
    %% Fourier-Bessel transform, get electrical field at z.
    % note: there is no diffrence between the transform and
    % inverse-transform operations. Refer to the book "Fourier Optics" by
    % Joseph W. Goodman, Third Edition, page 10, Chapter 2.1.5
    maskZ=zeros(size(r));
    for ii=1:length(r)
        a=Fz.*besselj(0,2*pi*r(ii)*rho)*2*pi.*rho;
        maskZ(ii)=trapz(rho,a);
    end
    trash=0;    
end

