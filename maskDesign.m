function [output]=maskDesign(beamD,f1,axicon,wavelength)
% MASK_DESIGN suggests 4 masks to be used and calculate the power loss of the mask; 
% the inner and the outer diameters of 4 masks are located at 1/e, 1/(5e), 
% 1/(10e) and 1/(15e) of the peak amplitude of the ring, respectively.
% the optical setup: axicon-->Lens 1 -->mask
% input argument beamD is the diameter of the Gaussian beam at axicon;
% f1 is the focal length of lens after the axicon;
%axicon includes the paraemters of axicon
%axicon.alpha is alpha angle of axicon which is 0.5*(180-apexAngle)
%axicon.refind is the refractive index of axicon material
% axicon.diameter is the diameter of the axicon
% options: https://www.asphericon.com/wp-content/uploads/SPA-North-America_2017.pdf
% the output contains the field:
%   output.thickness, thickness of the ring
%   output.outerD, outer diameter of the ring
%   output.innerD, inner diameter of the ring
%   output.field_mask, the eletrical field at the mask
%   output.r_mask, the coordinate at the mask
%   output.axicon, axicon information, save as the input axicon
%   output.wavelength, wavelength used, save as input wavelength
%   output.f1, focal length of the lens 1, save as input f1
%   output.beamD, beam diameter at the axicon, same as input beamD
%-------------example------------------
%     beamD=2.8; %the diameter of the Gaussian beam at axicon
%     f1=80;% f1 is the focal length of lens after the axicon
%     axicon.alpha=1; % alpha angle of axicon;
%     axicon.refind=1.4512;% refractive index of axicon material at the used wavelength; at 0.94
%     axicon.diameter=25;% mm
%     wavelength=0.94;
%    [output]=maskDesign(beamD,f1,axicon,wavelength)
if nargin==0
    %%
    beamD=2.8;%the diameter of the Gaussian beam at axicon
    f1=80;% f1 is the focal length of lens after the axicon
    axicon.alpha=1; % alpha angle of axicon;
    axicon.refind=1.4512;% refractive index of axicon material at the used wavelength; at 0.94
%     axicon.refind=1.4508; % at 970nm
    axicon.diameter=25;% mm
    wavelength=0.94;
% options: https://www.asphericon.com/wp-content/uploads/SPA-North-America_2017.pdf
%: 0.5, 1, 2, 5, 10, 20
end

%% simulation parameters at the plane of the axicon
step_axicon=0.4;% um, step size of each pixel for the coordinate on the axicon
pixelNumber=axicon.diameter*1000/step_axicon;% pixel number for the coordinate on the axicon
halfPixelNumber=pixelNumber*1/2;% pxiel numbers from center of the axicon to the edge
beta=(axicon.refind-1)*axicon.alpha/180*pi;% deflected angle of the beam
middleD=2*beta*f1;

%% simulation parameters at mask plane
step_mask=0.4;% micron;
r_mask=0:step_mask:middleD*1000*1; % radial coordinate at mask
initialPhase=0; % 
beamr=beamD/2*1000; % obtain beam radius at the axicon in the unit of um
r_axicon=(0:1:halfPixelNumber-1)*step_axicon;% radial coordinate at axicon

%% field at the axicon
% note: for the convex axicon there is a minus sign in the equation. we can
% refer to the function of the lens analogously which is
% exp(-1i*2*pi*r^2/2/f), where f is the focal length of the lens. 
field_axicon=exp(-1i*(2*pi*r_axicon*beta/wavelength-initialPhase));

%% Gaussian beam
Gauss=exp(-r_axicon.^2/beamr/beamr);
field_axicon=field_axicon.*Gauss;% phase added upon the Gaussian beam

%% field at the mask
field_mask=fourier_circular(field_axicon,wavelength,r_axicon,f1,r_mask);
amplitude=abs(field_mask);
ampRatio=1./([1,5,10,15]*exp(1));
[~,halfThickness]=findXatAmpRatio(amplitude,ampRatio);
halfThickness=halfThickness*step_mask/1000; % turn the unit into mm
outerD=2*(middleD/2+halfThickness);
innerD=2*(middleD/2-halfThickness);
residualPower_mask=zeros(size(outerD));
result='maskInfo';
codePath=mfilename('fullpath');
filePath=fileparts(codePath);
% filePath=cd;
resultPath=fullfile(filePath,result);
if ~exist(resultPath,'dir')
    mkdir(resultPath)
end

for ii=1:length(outerD)
    figure(ii);clf;
     residualPower_mask(ii)=mask_chracterization_plot(r_mask,field_mask,innerD(ii),outerD(ii));
     print('-r200',fullfile(resultPath,['f1_',num2str(f1),'mm_mask_',num2str(ii),'.png']),'-dpng')  
end

output.thickness=2*halfThickness;
output.outerD=outerD;
output.innerD=innerD;
output.field_mask=field_mask;
output.r_mask=r_mask;
output.axicon=axicon;
output.wavelength=wavelength;
output.f1=f1;
output.beamD=beamD;
save(fullfile(resultPath,['maskOutput_f1_',num2str(f1),'mm.mat']));

function residualPower_mask=mask_chracterization_plot(r_mask,field_mask,innerD,outerD)


%% calculate Transmittance after the mask
powerTmp=abs(field_mask).^2.*r_mask;
% flag1=r_mask*mp<=rOB*1000;%coordinates that are inside the objective  
flag2=(r_mask>=innerD/2*1000)&(r_mask<=outerD/2*1000); % coordinates that are inside the annular filter    
residualPower_mask=sum(powerTmp(flag2))/sum(powerTmp);
residualPower_mask=residualPower_mask*100;
%%
amplitude=abs(field_mask);
amplitude=amplitude/max(amplitude);
panelM=2;
panelN=1;
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
hAmp=plot(r_mask/1000,amplitude,'b','LineWidth',lineW);

limab=[0,0.5*(innerD+outerD)];

r2=outerD/2;
r1=innerD/2;
t1=text(0.1,0.5,['innerD=',num2str(r1*2,'%.3f'),'mm']);
t2=text(0.1,0.3,['outerD=',num2str(r2*2,'%.3f'),'mm']);


xlim(limab);
xlabel(['Transmittance=',num2str(residualPower_mask,'%.1f'),'%'])

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
angleY=angle(field_mask)/pi*180;
angleY(abs(angleY-360+180)<=10^-6)=-180;
ylabel('Amplitude','Fontname','Arial','Fontsize',labelSize)
xlim(limab);
subplot(panelM,panelN,plotOrder(plot_i));
plot_i=plot_i+1;
plot(r_mask/1000,angleY,'b','LineWidth',lineW);
% set(gca,'xTick',[0:0.3:1.8])
ylim([-200,200])
b=200;
hold on; 

%%
if r1>0.01
    plot([r1,r1],[-1,1]*b,'k--','LineWidth',lineW-0.5);
    hLine1=plot([r2,r2],[-1,1]*b,'k--','LineWidth',lineW-0.5);
    hLegend=legend([hLine1],'Annular mask');
    hLegend.Position=[0.1471 0.7964 0.1990 0.0501];    

end

box off;
set(gca,'FontName','Arial')
set(gca,'color','none')
set(gca,'TickDir','out');
set(gcf,'PaperPositionMode','auto')
set(gca,'yTick',[-180,0,180])
ylabel('Phase(Degree)','Fontname','Arial','Fontsize',labelSize)
xlim(limab);
xlabel('r(mm)','Fontname','Arial','Fontsize',labelSize);

function field_mask=fourier_circular(mask,wavelength,r_SLM,f,r_mask)
f=f*1000;
p=length(r_mask);  
field_mask=zeros(size(r_mask));
for ii=1:p
    a=mask.*besselj(0,2*pi*r_SLM*r_mask(ii)/wavelength/f)*2*pi.*r_SLM;
    field_mask(ii)=trapz(r_SLM,a);
end

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