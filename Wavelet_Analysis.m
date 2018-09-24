%  IMPORTANT!!!
%  remember to choose correct sum over divergence:  on-axis only, linear
%  sum (i.e. slit filter) or full radial sum


clear
path0 = './';

au_to_eV = 27.212;
eV_to_au = 1/au_to_eV;

W = importdata([path0,'frequency_p.dat']);      % (positive) angular frequency in eV
w = [W(1:end) ; -W(end-1:-1:2)] * eV_to_au;     % full angular fequency axis in au.
clear('W')

r = importdata([path0,'radius.dat']);
dr = r(2)-r(1);

P  = importdata([path0,'far_Propagated.dat']);
Parg  = importdata([path0,'far_argPropagated.dat']);

RN = 1;
% Hf_onaxis = sqrt([P(:,RN) ; flipud(P(2:(end-1),RN)) ]) .* [exp(1i*pi*Parg(:,RN)) ; exp(-1i*pi*flipud(Parg(2:(end-1),RN))) ] ;

sP = size(P);
Hf_linsum = zeros(2*(sP(1)-1),1);
for j = 1:length(Hf_linsum)
    if j <= sP(1)
        Hf_linsum(j) = sum(sqrt(P(j,:)).*exp(1i*pi*Parg(j,:)));
    else
        Hf_linsum(j) = sum(sqrt(P(length(Hf_linsum)+1-j,:)).*exp(-1i*pi*Parg(length(Hf_linsum)+1-j,:)));
    end
end
Hf_linsum = Hf_linsum*dr;

% Hf_radsum = zeros(2*(sP(1)-1),1);
% for j = 1:length(Hf_radsum)
%     if j <= sP(1)
%         Hf_radsum(j) = dr/16*sqrt(P(j,1))*exp(1i*pi*Parg(j,1)) + sum(sqrt(P(j,2:end)).*exp(1i*pi*Parg(j,2:end)).*r(2:end)');
%     else
%         Hf_radsum(j) = dr/16*sqrt(P(length(Hf_radsum)+1-j,1))*exp(-1i*pi*Parg(length(Hf_radsum)+1-j,1)) + sum(sqrt(P(length(Hf_radsum)+1-j,2:end)).*exp(-1i*pi*Parg(length(Hf_radsum)+1-j,2:end)).*r(2:end)');
%     end
% end
% Hf_radsum = Hf_radsum*dr*2*pi;
clear('P','Parg')

% Hf = Hf_onaxis;
Hf = Hf_linsum;
% Hf = Hf_radsum;


Nt = length(Hf);

T = 2*pi/(w(2)-w(1));
dt = T/Nt;
t = (dt-T/2) : dt : T/2;

tS = floor(length(t)*0.3);
tF = ceil(length(t)*0.7);
tJ = 4;
t2 = t(tS:tJ:tF);

max_eV=180;                                 % define max & min harmonics you want to look at
min_eV=15;

highW = max_eV*eV_to_au;
lowW = min_eV*eV_to_au;

%% test function
%Hf = fft(cos(2*t)+cos(3*t))';
%%

maxW = highW/2;     % max frequency
maxT = 2/Nt;        % max period                       %2*dt/(lowW*(t(end)-t(1)));

Na = 300;
sigma = 15%25;         % wavelet sigma parameter - central frequency = sigma / a

nS=log2(sigma/(2.0*maxW));
nF=log2(maxT*Nt*sigma/(2.0));

dn=(nF-nS)/Na;

count = 0;
a = zeros(1,Na);
integral = zeros(Na,Nt);

for n = nS:dn:nF
    count = count+1;
    a(count) = 2^n;

    wavelet = sqrt(a(count))*pi^-0.25 * exp(-(sigma - a(count)*w).^2/2);
    
    integral(count,:) = sqrt(2*pi) * ifft(Hf.*wavelet);
end

omega = (sigma./a)*au_to_eV;

figure; clf
% pcolor(0.024189/2.5*t,omega,abs(integral).^2);
% ylabel('Harmonic Energy / eV ');
pcolor(0.024189*t2,omega,abs(integral(:,tS:tJ:tF)).^2); shading flat
ylabel('Harmonic Energy / eV');

clmx = max(max(abs(integral(omega>50,:)).^2));
caxis([clmx/50 clmx])

xlim([-3 3])
xlabel('Time [cycles]');

%%
% 
% %CM = colormap(flipud(hot(256)));
% CMn = 256;
% jw1 = 0.2;
% jw2 = 0.95;
% JW = ones(CMn , 3);
% 
% tmp = floor(CMn*jw1):floor(CMn*jw2);
% JW(tmp , 1) = 2/3 - (2/3)*( tmp   - tmp(1)) / (tmp(end)  - tmp(1));
% JW(tmp, 2 ) = 1 - 0.4*exp( -(  10*(tmp-tmp(end - floor(length(tmp)/2)))/tmp(end) ).^2  );
% 
% tmp = 1 : floor(CMn*jw1);
% JW(tmp , 1) = 2/3;
% JW(tmp , 2) = (tmp - 1)/ tmp(end) ;
% 
% tmp = floor(CMn*jw2):CMn;
% JW(tmp, 3 ) = 1 - 0.6*(tmp - tmp(1))/(tmp(end)-tmp(1));
% 
% JW = hsv2rgb(JW);
% 
% colormap(JW);
% %colormap(flipud(hot(CMn)));

%%

c=0:.001:1;

rainbow_white(:,1) = max(min(1,23/9-50/9*c),max(0,50/9*c-41/9));
rainbow_white(:,2) = max(1-10*c,min(min(50/9*c-5/9,1),max(41/9-50/9*c,0)));
rainbow_white(:,3) = max(1-10*c,min(max(0,50/9*c-23/9),1));

colormap(rainbow_white)

set(gca,'xlim',[-10.9 10.9],'xtick',-10:2:10,'ylim',[26.8 129],'ytick',30:20:110)
% set(gca,'clim',[6e-26 8e-23])
set(gca,'clim',[3e-36 3.e-33])

%%


CB = colorbar;
set(get(CB,'children'),'ydata',[0 1])
%set(CB,'yscale','log')
%set(CB,'yticklabel',[0 0.25 0.5 0.75 1])
set(CB,'ylim',[0 1])
ylabel('Normalised magnitude of wavelet transform','parent',CB)
