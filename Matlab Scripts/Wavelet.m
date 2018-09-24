%  IMPORTANT!!!
%  remember to choose correct sum over divergence:  on-axis only, linear
%  sum (i.e. slit filter) or full radial sum


% clear

path=strcat('../');

spec = importdata([path,'Spec.dat']);  
omega0 = spec(1,2);

au_to_eV = 27.212;
eV_to_au = 1/au_to_eV;

clear('w'); w = importdata([path,'AxisW.dat']);     
% clear('r'); r = importdata([path,'AxisR.dat']);
% dr = r(2)-r(1);
clear('r'); r = importdata([path,'AxisRp.dat']);
dr = r(2,1)-r(1,1);
ddiv = r(2,2)-r(1,2);

clear('PX','PargX','HX');
PX = importdata([path,'HarmFarRWnormX.dat']);
PargX = -importdata([path,'HarmFarRWphaseX.dat']);
% PX = importdata([path,'HarmNearRWnormX.dat']);
% PargX = -importdata([path,'HarmNearRWphaseX.dat']);

opHarmInt = 'o'
% on-axis
if opHarmInt == 'o'
    RN = 1;
    RN = ceil(8e-3/ddiv);   % cut at 4mrad
%     HfX = sqrt(PX(:,RN)) .* exp(1i*PargX(:,RN));
    HfX = sqrt(PX(1:end,RN)) .* exp(1i*PargX(1:end,RN));
end
% linear sum
if opHarmInt == 'l'
    sP = size(PX);  % PY should always be the same size as PX
    HfX = zeros(sP(1),1);
    for j = 1:length(HfX)
        HfX(j) = sum(sqrt(PX(j,:)).*exp(1i*PargX(j,:)));
    end
    HfX = HfX*dr;
end
% radial integration
if opHarmInt == 'r'
    sP = size(PX);
    HfX = zeros(sPX(1),1);
    for j = 1:length(HfX)
        HfX(j) = dr/16*sqrt(PX(j,1))*exp(1i*PargX(j,1)) + sum(sqrt(PX(j,2:end)).*exp(1i*PargX(j,2:end)).*r(2:end)');
    end
    HfX = HfX*dr*2*pi;
end
clear('PX','PargX')





w = w(1:end/2);
HfX = HfX(1:end/2,:);





Nt = length(HfX);

T = 2*pi/(w(2)-w(1));
dt = T/Nt;
t = (dt-T/2) : dt : T/2;

tS = floor(length(t)*0.01);
tF = ceil(length(t)*0.99);
tJ = 1%4;
t2 = t(tS:tJ:tF);

max_eV=70;                                 % define max & min harmonics you want to look at
min_eV=20;

highW = max_eV*eV_to_au;
lowW = min_eV*eV_to_au;

maxW = highW/2;     % max frequency
maxT = 2/Nt/lowW;        % max period                       %2*dt/(lowW*(t(end)-t(1)));

Na = 300;
sigma = 13%25;      % wavelet sigma parameter - central frequency = sigma / a

nS=log2(sigma/(2.0*maxW));
nF=log2(maxT*Nt*sigma/(2.0));

dn=(nF-nS)/Na;
% dn=(nF-nS)/(Na-1);

count = 0;
a = zeros(1,Na);
clear('integralX')
integralX = zeros(Na,Nt);

for n = nS:dn:nF
    count = count+1;
    a(count) = 2^n;

%     Luke's formula - doesn't seem to conserve area under wavelet
%     pref = 4*sqrt(pi)/(1+erf(sigma))*a(count);
%     wavelet = pref*exp(-(sigma - a(count)*w).^2/2);

    wavelet = sqrt(a(count))*pi^-0.25 * exp(-(sigma - a(count)*w).^2/2);
    
%     integralX(count,:) = sqrt(2*pi) * ifft(HfX.*wavelet);
    
    integralX(count,:) = sqrt(2*pi) * fftshift(ifft(HfX.*wavelet));
    
    
end
clear('HfX')

omega = (sigma./a)*au_to_eV;

%%
% clmxX = max(max(abs(integralX(omega>1,:)).^2));
clmxX = max(max(abs(integralX(omega>20,tS:tJ:tF)).^2));

% figure; clf
% % pcolor(0.024189*t2,omega,abs(integral(:,tS:tJ:tF)).^2); shading flat; xlabel('Time / fs','fontsize',18)
% pcolor(t2*omega0/(2*pi),omega,abs(integral(:,tS:tJ:tF)).^2/clmx); shading flat; xlabel('Time / laser cycles','fontsize',18)
% ylabel('Photon Energy / eV','fontsize',18);

scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4)*1/2 scrsz(3)*.5 scrsz(4)*1/2])
pcolor(t2*omega0/(2*pi),omega,abs(integralX(:,tS:tJ:tF)).^2/clmxX); shading flat; 
xlabel('Time / laser cycles','fontsize',18);  ylabel('Photon Energy / eV','fontsize',18);
% set(gca,'xlim',[-1.1 3.1],'xtick',-3:.5:3,'ylim',[20 60],'ytick',20:10:200,'clim',[0 1],'fontsize',18)
CB = colorbar;
ylabel('Normalised Intensity','parent',CB,'fontsize',18)

% c=0:.001:1;
% rainbow_white(:,1) = max(min(1,23/9-50/9*c),max(0,50/9*c-41/9));
% rainbow_white(:,2) = max(1-10*c,min(min(50/9*c-5/9,1),max(41/9-50/9*c,0)));
% rainbow_white(:,3) = max(1-10*c,min(max(0,50/9*c-23/9),1));
% colormap(rainbow_white)


%% attopulse reconstruction

% P = importdata([path,'HarmFarRWnormX.dat']);
% Parg = importdata([path,'HarmFarRWphaseX.dat']);
% Hf = sqrt(P(:,RN)) .* exp(1i*Parg(:,RN));
% Hf_filt = Hf; Hf_filt(abs(w*27.212)<35) = 0; Hf_filt(abs(w*27.212)>60) = 0;
% DA = ifft(Hf);
% DA_filt = ifft(Hf_filt);
% figure; semilogy(w*27.212,abs(Hf(:,1)).^2)
% hold on; semilogy(w*27.212,abs(Hf_filt(:,1)).^2,'color','r')
% figure; plot(t*0.024189,real(DA))
% hold on; plot(t*0.024189,real(DA_filt),'color','r')

