%% notes:
% thin on-axis interference in far-field -> increase tRange

%%  control options
opPulseEtX = true;
opPulseEwX = false%true;
opPulseErzX = false%true;
opPulseErwX = false%true;
opPulseErzwX = false;
opPulseErtX = false;
opPulseErztX = true% false;
opPulseEztX = false%true;

opHarmNearErwX = true;
opHarmNearErzwX = false%true;
opHarmFarErwX = true;

opPulseEtY = false%true;
opPulseEwY = false;
opPulseErzY = false;
opPulseErwY = false;
opPulseErzwY = false;
opPulseErtY = false%true;
opPulseErztY = false;
opPulseEztY = false;

opHarmNearErwY = true;
opHarmNearErzwY = false;
opHarmFarErwY = true;

opHarmFarErwXY = false%true;  maxExy = 2.65e-19

%%  initial governing and axes variables
scrsz = get(0,'ScreenSize');

path='../';

spec = importdata([path,'Spec.dat']);
Nr = spec(2,4);
Nz = spec(2,5);
Ntp = spec(2,6);

au_to_m = 5.29177e-11;
m_to_au = 1.0/au_to_m;

mpi_rank = spec(3,1);
mpi_size = spec(3,2);

Dim = spec(1,6);

clear('r','z','tp','w','r2','z2','tp2','w2')
r = importdata([path,'AxisR.dat']);
z = importdata([path,'AxisZ.dat']);
tp = importdata([path,'AxisTp.dat']);
w = importdata([path,'AxisW.dat']);
r2 = [-flipud(r); r]*au_to_m*1e6;
z2 = z*au_to_m*1e3;
tp2 = tp*2.4189e-17*1e15;
w2 = w*27.212;
clear('r','z','tp','w')

rp = importdata([path,'AxisRp.dat']);
div = [-flipud(rp(:,2)); rp(:,2)]*1e3;
clear('rp')

%%  X-COMPONENTS

%%  plots on-axis laser field polarised along x-axis at the focus vs time
if opPulseEtX
    clear('Etp','Etp2')
    Etp = importdata([path,'PulseTpX.dat']);
    Etp2 = Etp*5.142e11*1e-9;
    clear('Etp')
    figure 
    plot(tp2, Etp2)
    xlabel('t'' / fs','fontsize',18); ylabel('E(t'') / Vnm^{-1}','fontsize',18)
end

%%  plots on-axis laser field polarised along x-axis in frequency domain at the focus vs frequency
if opPulseEwX
    clear('Ew')
    Ew = importdata([path,'PulseWnormX.dat']);
    figure 
    plot(w2, Ew)
    xlabel('\omega / eV','fontsize',18); ylabel('E(w)','fontsize',18); 
end

%%  plots spatial laser field polarised along x-axis at t0
if opPulseErzX
    clear('Erz')
    Erz = importdata([path,'PulseRZX.dat']);
    maxE = max(max(abs(Erz)))   % maxE = max(max(Erz))
    Erz2 = [flipud(Erz'); Erz'];
    clear('Erz')
    figure 
    imagesc(z2, r2, (Erz2/maxE).^2); set(gca,'ydir','normal','fontsize',18)
    set(gca,'ylim',[-100 100])
    xlabel('z / mm','fontsize',18); ylabel('r / {\mu}m','fontsize',18)
end

%%  plots laser field polarised along x-axis vs radius and frequency after interaction region
if opPulseErwX
    clear('E','E2')
    E = importdata([path,'PulseRWnormX.dat']);
    maxE = max(max(E));
    E2 = [flipud(E'); E'];
    clear('E')
    % plots (w,r,E) running through z
    figure('Position',[1 scrsz(4)*1/2 scrsz(3)*.5 scrsz(4)*1/2])
    imagesc(w2(1:Ntp/2), r2, (E2(:,1:Ntp/2)/maxE).^2); set(gca,'ydir','normal','fontsize',18); 
    xlabel('\omega / eV','fontsize',18); ylabel('r / {\mu}m','fontsize',18)
    set(gca,'clim',[0 1]); colorbar
    % figure('Position',[1 scrsz(4)*1/2 scrsz(3)*.5 scrsz(4)*1/2])
    % imagesc(w2(1+Ntp/2:Ntp), r2, (E2(:,1+Ntp/2:Ntp)/maxE).^2); set(gca,'ydir','normal','fontsize',18); 
    % xlabel('\omega / eV','fontsize',18); ylabel('r / {\mu}m','fontsize',18)
    % set(gca,'clim',[0 1]); colorbar
end

%%  plots laser field polarised along x-axis vs radius and frequency through interaction region as an animation
if opPulseErzwX
    clear('F','E_frag','E','E2')
    F(Ntp) = struct('cdata',[],'colormap',[]);

    % creates a cell array (signified by {}) with each cell containing 2d double array (tp,z vs r)
    for j = 1:mpi_size
        E_frag{j} = importdata(strcat(path,'PulseRZWnormX_C',num2str(j-1),'.dat'));
    end

    size_E_frag = size(E_frag{1});

    filter = 4;
    for j = 1:Nz
        for k=1:mpi_size
            E(1+(j-1)*(Ntp/filter)+(k-1)/mpi_size*(Ntp/filter):(j-1)*(Ntp/filter)+k/mpi_size*(Ntp/filter),1:size_E_frag(2)) = E_frag{k}(1+(j-1)*(Ntp/filter)/mpi_size:j*(Ntp/filter)/mpi_size,1:size_E_frag(2));
        end
    end
    clear('E_frag')

    maxE = max(max(E));
    E2 = [flipud(E'); E'];
    clear('E')

    mkdir video

    % plots (w,r,E) running through z for only positive w
    maxX=3e-4;
    counter = 0;
    filter = 4;
    figure('Position',[1 scrsz(4)*1/2 scrsz(3)*.5 scrsz(4)*1/2])
    for j = 1:4:Nz 
        counter = counter+1;
        %     full field
        imagesc(w2(1:filter:Ntp/2), r2, (E2(:,1+(j-1)*(Ntp/filter):j*(Ntp/filter)-(Ntp/filter)/2)/maxE).^2); set(gca,'ydir','normal','fontsize',18); 
        %     just one E_frag
        %     imagesc(w2(1:Ntp/mpi_size/2), r2, (E2(:,1+(j-1)*Ntp/mpi_size:j*Ntp/mpi_size-Ntp/mpi_size/2)/maxE).^2); set(gca,'ydir','normal','fontsize',18); 
        xlabel('\omega / eV','fontsize',18); ylabel('r / {\mu}m','fontsize',18)
        %     set(gca,'xlim',[0 maxX]); 
        set(gca,'clim',[0 1]); colorbar
        text(maxX*0.7,max(r2)*0.8,['z = ' num2str(z2(j),'%4.1f') 'mm'],'Color','w','FontSize',18)
        F(counter) = getframe(gcf);
        imwrite(F(counter).cdata, sprintf('./video/%04d.png',counter));
    end
end

%%  plots laser field polarised along x-axis vs radius and time after interaction region
if opPulseErtX
    clear('E','E2')
    E = importdata([path,'PulseRTpX.dat']);
    maxE = max(max(E));
    E2 = [fliplr(E') E'];
    clear('E')
    % plots (t',r,E) running through z
    figure('Position',[1 scrsz(4)*1/2 scrsz(3)*.5 scrsz(4)*1/2])
    imagesc(tp2, r2, (E2'/maxE).^2); set(gca,'ydir','normal','fontsize',18)
    xlabel('t'' / fs','fontsize',18); ylabel('r / {\mu}m','fontsize',18)
    set(gca,'clim',[0 1]); colorbar
end

%%  plots laser field polarised along x-axis vs radius and time through interaction region as an animation
if opPulseErztX
    clear('F','E_frag','E','E2')
    F(Ntp) = struct('cdata',[],'colormap',[]);
    % creates a cell array (signified by {}) with each cell containing 2d double array (tp,z vs r)
    for j = 1:mpi_size
        E_frag{j} = importdata(strcat(path,'PulseRZTpX_C',num2str(j-1),'.dat'));
    end
    size_E_frag = size(E_frag{1});

    for j = 1:Nz
        for k=1:mpi_size
            E(1+(k-1)*(Nr/mpi_size):k*(Nr/mpi_size),1:size_E_frag(2)) = E_frag{k}(1+(j-1)*(Nr/mpi_size):j*(Nr/mpi_size),1:size_E_frag(2));
        end
        E2(1+(j-1)*(2*Nr):j*(2*Nr),1:size_E_frag(2)) = [flipud(E); E];
    end
    clear('E_frag')

    maxE = max(max(E2));
    clear('E')

    mkdir video

    % plots (t',r,E) running through z
    filter = spec(3,5);%4;
    figure('Position',[1 scrsz(4)*1/2 scrsz(3)*.5 scrsz(4)*1/2])
    frame=0;
    tlim=10;
    rlim=100;
    for j = 1:4:Nz 
        frame=frame+1;
        imagesc(tp2(1:filter:end), r2, (E2(1+(j-1)*(2*Nr):j*(2*Nr),:)/maxE).^2); set(gca,'ydir','normal','fontsize',18)
        xlabel('t'' / fs','fontsize',18); ylabel('r / {\mu}m','fontsize',18)
        set(gca,'xlim',[-tlim tlim]);
        set(gca,'ylim',[-rlim rlim]);
        set(gca,'clim',[0 1]); colorbar
        text(tlim*0.6,rlim*0.8,['z = ' num2str(z2(j),'%4.1f') 'mm'],'Color','w','FontSize',18)
        F(frame) = getframe(gcf);
        imwrite(F(frame).cdata, sprintf('./video/%04d.png',frame));
    end

    % plots on-axis (t',E) at the centre of gas jet (need to adjust j value if not at centre of z range)
    figure('Position',[1 scrsz(4)*1/2 scrsz(3)*.5 scrsz(4)*1/2])
    for j = 1+Nz/2 
        plot(tp2(1:filter:end), (E2(1+(j-1)*(2*Nr)+Nr,:)/maxE).^2); set(gca,'ydir','normal','fontsize',18)
        xlabel('t'' / fs','fontsize',18); ylabel('E','fontsize',18)
    end

    % % plots (z,r,E) running through t'
    % figure('Position',[1 scrsz(4)*1/2 scrsz(3)*.5 scrsz(4)*1/2])
    % frame=0;
    % for j = 1+Ntp/2  %Ntp/4:4:Ntp  %Ntp/16*7:1:Ntp/16*9 
    %     frame=frame+1;
    %     E = reshape(E2(:,j),2*Nr,[]);
    %     imagesc(z2, r2, (E/maxE).^2); set(gca,'ydir','normal','fontsize',18)
    %     xlabel('z / mm','fontsize',18); ylabel('r / {\mu}m','fontsize',18)
    %     set(gca,'clim',[0 1]); colorbar
    %     text(max(tp2)*0.7,max(r2)*0.8,['t'' = ' num2str(tp2(j),'%4.1f') 'fs'],'Color','w','FontSize',18)
    %     F(frame) = getframe(gcf);
    %     imwrite(F(frame).cdata, sprintf('./video/%04d.png',j));
    % end
end

%% plot on-axis laser field polarised along x-axis vs time through interaction region 
if opPulseEztX
    clear('E2')
    E2 = importdata([path,'PulseZTpX.dat']);
    maxE = max(max(E2));
%     % plots (t',r,E) running through z
%     figure('Position',[1 scrsz(4)*1/2 scrsz(3)*.5 scrsz(4)*1/2])
%     imagesc(tp2, z2, (E2/maxE).^2); set(gca,'ydir','normal','fontsize',18)
%     xlabel('t'' / fs','fontsize',18); ylabel('z / mm','fontsize',18)
%     set(gca,'clim',[0 1]); colorbar
    % animation version
    clear('F')
    F(Nz) = struct('cdata',[],'colormap',[]);
    mkdir video
    % plots (t',E) running through z
    filter = 1;
    figure('Position',[1 scrsz(4)*1/2 scrsz(3)*.5 scrsz(4)*1/2])
    frame=0;
    tlim=10;
    for j = 1:Nz 
        frame=frame+1;
        plot(tp2(1:filter:end), 5.142e2*E2(j,1:filter:end),'linewidth',2)
        xlabel('t'' / fs','fontsize',18); ylabel('E / Vnm^{-1}','fontsize',18)
        set(gca,'xlim',[-20 20],'fontsize',18);
        set(gca,'ylim',[-10*ceil(51.42*maxE) 10*ceil(51.42*maxE)]);
        text(tlim*0.7,8*ceil(51.42*maxE),['z = ' num2str(z2(j),'%4.1f') 'mm'],'Color','black','FontSize',18)
        F(frame) = getframe(gcf);
        imwrite(F(frame).cdata, sprintf('./video/%04d.png',frame));
    end
end

%%  plots harmonic field polarised along x-axis vs radius and frequency after interaction region
if opHarmNearErwX
    clear('E','E2')
    E = importdata([path,'HarmNearRWnormX.dat']);
    % maxE = max(max(E));
    maxE = max(max(E(abs(w2)>20,:)));
    E2 = [flipud(E'); E'];
    clear('E')
    % plots (w,r,E) running through z
    fltr = 8;
    figure('Position',[1 scrsz(4)*1/2 scrsz(3)*.5 scrsz(4)*1/2])
    imagesc(w2(1:fltr:Ntp/2), r2, E2(:,1:fltr:Ntp/2)/maxE); set(gca,'ydir','normal','fontsize',18); 
    xlabel('\omega / eV','fontsize',18); ylabel('r / {\mu}m','fontsize',18)
    set(gca,'xlim',[15 75])
    set(gca,'ylim',[-30 30])
    set(gca,'clim',[0 1]); colorbar
    
%     figure('Position',[1 scrsz(4)*1/2 scrsz(3)*.5 scrsz(4)*1/2])
%     semilogy(w2(1:Ntp/2), E2(end/2+1:end/2+7,1:Ntp/2)/maxE);
%     xlabel('\omega / eV','fontsize',18);
%     set(gca,'xlim',[20 60],'ylim',[1e-2 2e0])
end
    
%%  plots harmonic field polarised along x-axis vs radius and time through interaction region as an animation
if opHarmNearErzwX
    clear('F','E_frag','E','E2')
    F(Ntp) = struct('cdata',[],'colormap',[]);

    % creates a cell array (signified by {}) with each cell containing 2d double array (tp,z vs r)
    for j = 1:mpi_size
        E_frag{j} = importdata(strcat(path,'HarmNearRZWnormX_C',num2str(j-1),'.dat'));
    end

    size_E_frag = size(E_frag{1});

    filter = spec(3,5);%4;
    for j = 1:Nz
        for k=1:mpi_size
            E(1+(j-1)*(Ntp/filter)+(k-1)/mpi_size*(Ntp/filter):(j-1)*(Ntp/filter)+k/mpi_size*(Ntp/filter),1:size_E_frag(2)) = E_frag{k}(1+(j-1)*(Ntp/filter)/mpi_size:j*(Ntp/filter)/mpi_size,1:size_E_frag(2));
        end
    end
    clear('E_frag')

    maxE = max(max(E));
    E2 = [flipud(E'); E'];
    clear('E')

    mkdir video

    % plots (w,r,E) running through z for only positive w
    maxX=3e-4;
    counter = 0;

    figure('Position',[1 scrsz(4)*1/2 scrsz(3)*.5 scrsz(4)*1/2])
    for j = 1:4:Nz 
        counter = counter+1;
    % %     full field
        imagesc(w2(1:filter:Ntp/2), r2, (E2(:,1+(j-1)*(Ntp/filter):j*(Ntp/filter)-(Ntp/filter)/2)/maxE).^2); set(gca,'ydir','normal','fontsize',18); 
    % %     just one E_frag
    %     imagesc(w2(1:Ntp/mpi_size/2), r2, (E2(:,1+(j-1)*Ntp/mpi_size:j*Ntp/mpi_size-Ntp/mpi_size/2)/maxE).^2); set(gca,'ydir','normal','fontsize',18); 
        xlabel('\omega / eV','fontsize',18); ylabel('r / {\mu}m','fontsize',18)
    %     set(gca,'xlim',[0 maxX]); 
        set(gca,'clim',[0 1]); colorbar
        text(maxX*0.7,max(r2)*0.8,['z = ' num2str(z2(j),'%4.1f') 'mm'],'Color','w','FontSize',18)
        F(counter) = getframe(gcf);
        imwrite(F(counter).cdata, sprintf('./video/%04d.png',counter));
    end
end

%%  plots harmonic field polarised along x-axis vs radius and frequency after farfield propagation
if opHarmFarErwX
    clear('E','E2')
    E = importdata([path,'HarmFarRWnormX.dat']);
    maxEx = max(max(E(abs(w2)>20,:)));
    Ex = [flipud(E'); E'];
    clear('E')
    % plots (w,r,E) running through z
    figure('Position',[1 scrsz(4)*1/2 scrsz(3)*.5 scrsz(4)*1/2])
    imagesc(w2(1:fltr:Ntp/2), div, Ex(:,1:fltr:Ntp/2)/maxEx); set(gca,'ydir','normal','fontsize',18); 
    xlabel('\omega / eV','fontsize',18); ylabel('Divergence / mrad','fontsize',18)
    set(gca,'xlim',[15 75])
    set(gca,'ylim',[-10 10])
    set(gca,'clim',[0 1]); colorbar
    
       
%     figure('Position',[1 scrsz(4)*1/2 scrsz(3)*.5 scrsz(4)*1/2])
%     imagesc(w2(1:Ntp/2), div, log10(Ex(:,1:Ntp/2)/maxEx)); set(gca,'ydir','normal','fontsize',18);
%     xlabel('\omega / eV','fontsize',18); ylabel('Divergence / mrad','fontsize',18)
%     set(gca,'xlim',[20 60])
%     set(gca,'ylim',[-10 10])
%     set(gca,'clim',[-3 0]); colorbar

% % % nwS = 1435;
% % % nwF = 1969%1610;
% % % rS = 148%161%148;
% % % rF = 160%
% % % 
% % % clear('int_r_E'); int_r_E = zeros(nwF-nwS+1,1);
% % % counter=0;
% % % ddiv = div(100)-div(99);
% % % for nw=nwS:nwF
% % %     counter=counter+1;
% % %     int_r_E(counter) = sum(Ex(rF:256,nw).*div(rF:256)) * 1e-3/5.29177e-11 * ddiv * 2*pi;
% % % end
% % % int_rw_E_long = sum(int_r_E)*(w2(2)-w2(1))/27.212
% % % 
% % % clear('int_r_E'); int_r_E = zeros(1610-1435+1,1);
% % % counter=0;
% % % ddiv = div(100)-div(99);
% % % for nw=nwS:nwF
% % %     counter=counter+1;
% % %     int_r_E(counter) = sum(Ex(128:rS,nw).*div(128:rS)) * 1e-3/5.29177e-11 * ddiv * 2*pi;
% % % end
% % % int_rw_E_short = sum(int_r_E)*(w2(2)-w2(1))/27.212

end

%%

%%  Y-COMPONENTS

%%  plots on-axis laser field polarised along y-axis at the focus vs time
if opPulseEtY
    clear('Etp','Etp2')
    Etp = importdata([path,'PulseTpY.dat']);
    Etp2 = Etp*5.142e11*1e-9;
    clear('Etp')
    figure 
    plot(tp2, Etp2)
    xlabel('t'' / fs','fontsize',18); ylabel('E(t'') / Vnm^{-1}','fontsize',18)
end

%%  plots on-axis laser field polarised along y-axis in frequency domain at the focus vs frequency
if opPulseEwY
    clear('Ew')
    Ew = importdata([path,'PulseWnormY.dat']);
    figure 
    plot(w2, Ew)
    xlabel('\omega / eV','fontsize',18); ylabel('E(w)','fontsize',18)
end

%%  plots spatial laser field polarised along y-axis at t0
if opPulseErzY
    clear('Erz')
    Erz = importdata([path,'PulseRZY.dat']);
    maxE = max(max(Erz))
    Erz2 = [flipud(Erz'); Erz'];
    clear('Erz')
    figure 
    imagesc(z2, r2, (Erz2/maxE).^2); set(gca,'ydir','normal','fontsize',18)
    set(gca,'ylim',[-100 100])
    xlabel('z / mm','fontsize',18); ylabel('r / {\mu}m','fontsize',18)
end

%%  plots laser field polarised along y-axis vs radius and frequency after interaction region
if opPulseErwY
    clear('E','E2')
    E = importdata([path,'PulseRWnormY.dat']);
    maxE = max(max(E));
    E2 = [flipud(E'); E'];
    clear('E')
    % plots (w,r,E) running through z
    figure('Position',[1 scrsz(4)*1/2 scrsz(3)*.5 scrsz(4)*1/2])
    imagesc(w2(1:Ntp/2), r2, (E2(:,1:Ntp/2)/maxE).^2); set(gca,'ydir','normal','fontsize',18); 
    xlabel('\omega / eV','fontsize',18); ylabel('r / {\mu}m','fontsize',18)
    set(gca,'clim',[0 1]); colorbar
    % figure('Position',[1 scrsz(4)*1/2 scrsz(3)*.5 scrsz(4)*1/2])
    % imagesc(w2(1+Ntp/2:Ntp), r2, (E2(:,1+Ntp/2:Ntp)/maxE).^2); set(gca,'ydir','normal','fontsize',18); 
    % xlabel('\omega / eV','fontsize',18); ylabel('r / {\mu}m','fontsize',18)
    % set(gca,'clim',[0 1]); colorbar
end

%%  plots laser field polarised along y-axis vs radius and frequency through interaction region as an animation
if opPulseErzwY
    clear('F','E_frag','E','E2')
    F(Ntp) = struct('cdata',[],'colormap',[]);

    % creates a cell array (signified by {}) with each cell containing 2d double array (tp,z vs r)
    for j = 1:mpi_size
        E_frag{j} = importdata(strcat(path,'PulseRZWnormY_C',num2str(j-1),'.dat'));
    end

    size_E_frag = size(E_frag{1});

    filter = 4;
    for j = 1:Nz
        for k=1:mpi_size
            E(1+(j-1)*(Ntp/filter)+(k-1)/mpi_size*(Ntp/filter):(j-1)*(Ntp/filter)+k/mpi_size*(Ntp/filter),1:size_E_frag(2)) = E_frag{k}(1+(j-1)*(Ntp/filter)/mpi_size:j*(Ntp/filter)/mpi_size,1:size_E_frag(2));
        end
    end
    clear('E_frag')

    maxE = max(max(E));
    E2 = [flipud(E'); E'];
    clear('E')

    mkdir video

    % plots (w,r,E) running through z for only positive w
    maxX=3e-4;
    counter = 0;
    filter = 4;
    figure('Position',[1 scrsz(4)*1/2 scrsz(3)*.5 scrsz(4)*1/2])
    for j = 1:4:Nz 
        counter = counter+1;
        %     full field
        imagesc(w2(1:filter:Ntp/2), r2, (E2(:,1+(j-1)*(Ntp/filter):j*(Ntp/filter)-(Ntp/filter)/2)/maxE).^2); set(gca,'ydir','normal','fontsize',18); 
        %     just one E_frag
        %     imagesc(w2(1:Ntp/mpi_size/2), r2, (E2(:,1+(j-1)*Ntp/mpi_size:j*Ntp/mpi_size-Ntp/mpi_size/2)/maxE).^2); set(gca,'ydir','normal','fontsize',18); 
        xlabel('\omega / eV','fontsize',18); ylabel('r / {\mu}m','fontsize',18)
        %     set(gca,'xlim',[0 maxX]); 
        set(gca,'clim',[0 1]); colorbar
        text(maxX*0.7,max(r2)*0.8,['z = ' num2str(z2(j),'%4.1f') 'mm'],'Color','w','FontSize',18)
        F(counter) = getframe(gcf);
        imwrite(F(counter).cdata, sprintf('./video/%04d.png',counter));
    end
end

%%  plots laser field polarised along y-axis vs radius and time after interaction region
if opPulseErtY
    clear('E','E2')
    E = importdata([path,'PulseRTpY.dat']);
    maxE = max(max(E));
    E2 = [fliplr(E') E'];
    clear('E')
    % plots (t',r,E) running through z
    figure('Position',[1 scrsz(4)*1/2 scrsz(3)*.5 scrsz(4)*1/2])
    imagesc(tp2, r2, (E2'/maxE).^2); set(gca,'ydir','normal','fontsize',18)
    xlabel('t'' / fs','fontsize',18); ylabel('r / {\mu}m','fontsize',18)
    set(gca,'clim',[0 1]); colorbar
end

%%  plots laser field polarised along y-axis vs radius and time through interaction region as an animation
if opPulseErztY
    clear('F','E_frag','E','E2')
    F(Ntp) = struct('cdata',[],'colormap',[]);
    % creates a cell array (signified by {}) with each cell containing 2d double array (tp,z vs r)
    for j = 1:mpi_size
        E_frag{j} = importdata(strcat(path,'PulseRZTpY_C',num2str(j-1),'.dat'));
    end
    size_E_frag = size(E_frag{1});

    for j = 1:Nz
        for k=1:mpi_size
            E(1+(k-1)*(Nr/mpi_size):k*(Nr/mpi_size),1:size_E_frag(2)) = E_frag{k}(1+(j-1)*(Nr/mpi_size):j*(Nr/mpi_size),1:size_E_frag(2));
        end
        E2(1+(j-1)*(2*Nr):j*(2*Nr),1:size_E_frag(2)) = [flipud(E); E];
    end
    clear('E_frag')

    maxE = max(max(E2));
    clear('E')

    mkdir video

    % plots (t',r,E) running through z
    filter = 4;
    figure('Position',[1 scrsz(4)*1/2 scrsz(3)*.5 scrsz(4)*1/2])
    frame=0;
    tlim=10;
    rlim=100;
    for j = 1:4:Nz 
        frame=frame+1;
        imagesc(tp2(1:filter:end), r2, (E2(1+(j-1)*(2*Nr):j*(2*Nr),:)/maxE).^2); set(gca,'ydir','normal','fontsize',18)
        xlabel('t'' / fs','fontsize',18); ylabel('r / {\mu}m','fontsize',18)
    %     set(gca,'xlim',[-tlim tlim]);
        set(gca,'ylim',[-rlim rlim]);
        set(gca,'clim',[0 1]); colorbar
        text(tlim*0.6,rlim*0.8,['z = ' num2str(z2(j),'%4.1f') 'mm'],'Color','w','FontSize',18)
        F(frame) = getframe(gcf);
        imwrite(F(frame).cdata, sprintf('./video/%04d.png',frame));
    end

    % plots on-axis (t',E) at the centre of gas jet (need to adjust j value if not at centre of z range)
    figure('Position',[1 scrsz(4)*1/2 scrsz(3)*.5 scrsz(4)*1/2])
    for j = 1+Nz/2 
        plot(tp2(1:filter:end), (E2(1+(j-1)*(2*Nr)+Nr,:)/maxE).^2); set(gca,'ydir','normal','fontsize',18)
        xlabel('t'' / fs','fontsize',18); ylabel('E','fontsize',18)
    end

    % % plots (z,r,E) running through t'
    % figure('Position',[1 scrsz(4)*1/2 scrsz(3)*.5 scrsz(4)*1/2])
    % frame=0;
    % for j = 1+Ntp/2  %Ntp/4:4:Ntp  %Ntp/16*7:1:Ntp/16*9 
    %     frame=frame+1;
    %     E = reshape(E2(:,j),2*Nr,[]);
    %     imagesc(z2, r2, (E/maxE).^2); set(gca,'ydir','normal','fontsize',18)
    %     xlabel('z / mm','fontsize',18); ylabel('r / {\mu}m','fontsize',18)
    %     set(gca,'clim',[0 1]); colorbar
    %     text(max(tp2)*0.7,max(r2)*0.8,['t'' = ' num2str(tp2(j),'%4.1f') 'fs'],'Color','w','FontSize',18)
    %     F(frame) = getframe(gcf);
    %     imwrite(F(frame).cdata, sprintf('./video/%04d.png',j));
    % end
end

%%  plots harmonic field polarised along y-axis vs radius and frequency after interaction region
if opHarmNearErwY
    clear('E','E2')
    E = importdata([path,'HarmNearRWnormY.dat']);
    % maxE = max(max(E));
    maxE = max(max(E(abs(w2)>20,:)));
    E2 = [flipud(E'); E'];
    clear('E')
    % plots (w,r,E) running through z
    figure('Position',[1 scrsz(4)*1/2 scrsz(3)*.5 scrsz(4)*1/2])
    imagesc(w2(1:Ntp/2), r2, E2(:,1:Ntp/2)/maxE); set(gca,'ydir','normal','fontsize',18); 
    xlabel('\omega / eV','fontsize',18); ylabel('r / {\mu}m','fontsize',18)
    set(gca,'xlim',[20 60])
    set(gca,'ylim',[-70 70])
    set(gca,'clim',[0 1]); colorbar
end
    
%%  plots harmonic field polarised along y-axis vs radius and time through interaction region as an animation
if opHarmNearErzwY
    clear('F','E_frag','E','E2')
    F(Ntp) = struct('cdata',[],'colormap',[]);

    % creates a cell array (signified by {}) with each cell containing 2d double array (tp,z vs r)
    for j = 1:mpi_size
        E_frag{j} = importdata(strcat(path,'HarmNearRZWnormY_C',num2str(j-1),'.dat'));
    end

    size_E_frag = size(E_frag{1});

    filter = 1;
    for j = 1:Nz
        for k=1:mpi_size
            E(1+(j-1)*(Ntp/filter)+(k-1)/mpi_size*(Ntp/filter):(j-1)*(Ntp/filter)+k/mpi_size*(Ntp/filter),1:size_E_frag(2)) = E_frag{k}(1+(j-1)*(Ntp/filter)/mpi_size:j*(Ntp/filter)/mpi_size,1:size_E_frag(2));
        end
    end
    clear('E_frag')

    maxE = max(max(E));
    E2 = [flipud(E'); E'];
    clear('E')

    mkdir video

    % plots (w,r,E) running through z for only positive w
    maxX=3e-4;
    counter = 0;

    figure('Position',[1 scrsz(4)*1/2 scrsz(3)*.5 scrsz(4)*1/2])
    for j = 1:4:Nz 
        counter = counter+1;
    % %     full field
        imagesc(w2(1:filter:Ntp/2), r2, (E2(:,1+(j-1)*(Ntp/filter):j*(Ntp/filter)-(Ntp/filter)/2)/maxE).^2); set(gca,'ydir','normal','fontsize',18); 
    % %     just one E_frag
    %     imagesc(w2(1:Ntp/mpi_size/2), r2, (E2(:,1+(j-1)*Ntp/mpi_size:j*Ntp/mpi_size-Ntp/mpi_size/2)/maxE).^2); set(gca,'ydir','normal','fontsize',18); 
        xlabel('\omega / eV','fontsize',18); ylabel('r / {\mu}m','fontsize',18)
    %     set(gca,'xlim',[0 maxX]); 
        set(gca,'clim',[0 1]); colorbar
        text(maxX*0.7,max(r2)*0.8,['z = ' num2str(z2(j),'%4.1f') 'mm'],'Color','w','FontSize',18)
        F(counter) = getframe(gcf);
        imwrite(F(counter).cdata, sprintf('./video/%04d.png',counter));
    end
end

%%  plots harmonic field polarised along y-axis vs radius and frequency after farfield propagation
if opHarmFarErwY
    clear('E','E2')
    E = importdata([path,'HarmFarRWnormY.dat']);
    % maxE = max(max(E));
    maxEy = max(max(E(abs(w2)>20,:)));
    Ey = [flipud(E'); E'];
    clear('E')
    % plots (w,r,E) running through z
    figure('Position',[1 scrsz(4)*1/2 scrsz(3)*.5 scrsz(4)*1/2])
    imagesc(w2(1:Ntp/2), div, Ey(:,1:Ntp/2)/maxEy); set(gca,'ydir','normal','fontsize',18); 
    xlabel('\omega / eV','fontsize',18); ylabel('Divergence / mrad','fontsize',18)
    set(gca,'xlim',[20 60])
    set(gca,'ylim',[-10 10])
    set(gca,'clim',[0 1]); colorbar
end

%%
if opHarmFarErwXY
    figure('Position',[1 scrsz(4)*1/2 scrsz(3)*.5 scrsz(4)*1/2])
    imagesc(w2(1:Ntp/2), div, (Ex(:,1:Ntp/2)+Ey(:,1:Ntp/2))/maxExy); set(gca,'ydir','normal','fontsize',18);
    xlabel('\omega / eV','fontsize',18); ylabel('Divergence / mrad','fontsize',18)
    set(gca,'xlim',[20 60])
    set(gca,'ylim',[-10 10])
    set(gca,'clim',[0 1]); colorbar
end

