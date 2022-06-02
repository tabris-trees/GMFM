clear;

%% Load the data
disp('Loading the data...');
datafolder = 'Raw_data\';
B_data = Read_data(datafolder);

% plot the position of stations
obs_theta = B_data(:,1);
obs_phi = B_data(:,2);
f = figure(5);
set(f,'Units','centimeters','Position',[10,5,25,18]);
clf;
m_proj('robinson','clongitude',180);
hold on;
m_coast('patch',[.9 .9 .9],'edgecolor','none'); % draw coastlines
s = m_plot(obs_phi,obs_theta,'^','markersize',5,'color','r','MarkerFaceColor','r');
m_grid('tickdir','in','backcolor',[.2 .65 1],'linewi',2); % 边框属性）
title('Observatories used');
legend(s,'Stations location','Location','southoutside');
print(f,'stations','-r600','-djpeg');


%% Inverse the Gauss Coefficient
disp('Computing the Gauss Coefficients...');
order = 4;
[Gauss_efficient,kernel_matrix] = Inverse(B_data,order);

%% Evaluate model error
disp('Evaluating the model error...')
est_obs_all = kernel_matrix*Gauss_efficient;
est_obs = reshape(est_obs_all,3,length(kernel_matrix)/3)';
obs_value = B_data(:,3:end);
f1 = vecnorm(est_obs,2,2);
f2 = vecnorm(obs_value,2,2);
% Draw a picture to illustrate the error
f = figure(16);
clf;
set(f,'Units','centimeters','Position',[10,5,25,18]);
errorf(1,f1,f2,'F');
errorf(2,est_obs(:,1),obs_value(:,1),'X');
errorf(3,est_obs(:,2),obs_value(:,2),'Y');
errorf(4,est_obs(:,3),obs_value(:,3),'Z');
print(f,'est_obs','-r600','-djpeg');
error = abs(f2-f1);
f = figure(17);
set(f,'Units','centimeters','Position',[10,5,25,18]);
clf;
m_proj('robinson','clongitude',180);
hold on;
m_coast('patch',[.9 .9 .9],'edgecolor','none'); % draw coastlines
s = m_scatter(obs_phi,obs_theta,error/10,'^','color','r','MarkerFaceColor','r');
m_grid('tickdir','in','backcolor',[.2 .65 1],'linewi',2); % 边框属性）
title('Errors of stations');
l = legend(s,'errors value','Location','southoutside');
print(f,'errors','-r600','-djpeg');

%% Compute the generated global model
disp('Computing the global model...')
% Gauss_efficient = importdata("Gauss_Coefficient.txt");order = 4;
est_theta = linspace(-90,90-0.05,81);
est_phi = linspace(-180,180,81);
Compute(est_theta,est_phi,Gauss_efficient,order);

%% Figure
disp('Figuring the global model...')
% est_theta = linspace(-90,90-0.05,81);est_phi = linspace(-180,180,81);
est_data = importdata("MMD_20200305_global.txt");

F = est_data(:,end); % F
Br = est_data(:,6); % Br
D = est_data(:,7)*(180/pi); % D
I = est_data(:,8)*(180/pi); % I

fig_value_c(est_theta,est_phi,F,11,'Total Field (F) in nT 20200305','F-c');
fig_value_c(est_theta,est_phi,Br,22,'Radial component (B_r) in nT 20200305','Br-c');
fig_value_c(est_theta,est_phi,D,33,'Declination (D) in degree 20200305','D-c');
fig_value_c(est_theta,est_phi,I,44,'Inclination (I) in degree 20200305','I-c');

fig_value_p(est_theta,est_phi,F,1,'F(nT)','F-p');
fig_value_p(est_theta,est_phi,Br,2,'B_r(nT)','Br-p');
fig_value_p(est_theta,est_phi,D,3,'D(degree)','D-p');
fig_value_p(est_theta,est_phi,I,4,'I(degree)','I-p');

%% comparing the spectra with that of IGRF-13
disp('Figuring the spectrum...')
% Gauss_efficient = importdata("Gauss_Coefficient.txt");order = 4;
GS_IGRF_all = importdata("igrf13coeffs.txt").data;
GS_IGRF = GS_IGRF_all(:,27);

spectrum_SH = Spectrum(Gauss_efficient,order);
spectrum_SH_IGRF = Spectrum(GS_IGRF,13);

f = figure(6);
set(f,'Units','centimeters','Position',[10,5,25,18]);
clf;
semilogy(spectrum_SH,'-o','LineWidth',1.5);
hold on;
semilogy(spectrum_SH_IGRF,'-p','LineWidth',1.5);
xlabel('Degree n');
ylabel('(nT^2)');
legend(['My ',num2str(order),' degree model'],'IGRF-13 model');
title('Comparing the spectra with that of IGRF-13');
print(f,'spectrum','-r600','-djpeg');

%% Functions
function fig_value_c(est_theta,est_phi,value_line,n,titles,name)
    [LG,LT] = meshgrid(est_phi,est_theta);
    value = reshape(value_line,length(est_phi),length(est_theta))';
    f = figure(n);
    set(f,'Units','centimeters','Position',[10,5,25,18]);
    clf;
    m_proj('Miller');
    hold on;
    % contour plot
    m_coast('patch',[.9 .9 .9],'LineWidth',1); % draw coastlines
    m_grid('tickdir','in','backcolor',[101,231,255]./255,'linewi',2); % 边框属性）
    [c0,h0] = m_contour(LG,LT,value,30,'LineWidth',1.5);
    clabel(c0,h0,'manual','FontWeight','bold','Color', ...
        'k','FontSize',10,'LabelSpacing',500); % Select the label manually
    mycmap = [0 0 1;1 0 0];
    if n == 11 || n == 22
        mycmap = [1,0,0]; % If draw a graph of F, use all the red lines
    end
    colormap(mycmap);
    colorbar('southoutside');
    title(titles);
    print(f,name,'-r600','-djpeg');
end

function fig_value_p(est_theta,est_phi,value_line,n,label,name)
    [LG,LT] = meshgrid(est_phi,est_theta);
    value = reshape(value_line,length(est_phi),length(est_theta))';
    f = figure(n);
    set(f,'Units','centimeters','Position',[10,5,25,18]);
    clf;
    m_proj('Hammer','clong',-360);
    % contourf plot
    hold on;
    m_pcolor(LG-360,LT,value);shading interp;
    m_coast('Color',[0 0 0],'LineWidth',1); % draw coastlines
    m_grid('xaxis','middle','linewi',2);
    colormap(getcolor('colormap.png'));
    h = colorbar('southoutside');
    set(get(h,'xlabel'),'string',label);
    print(f,name,'-r600','-djpeg');
end

function Sn=Spectrum(GS,order)
    Sn = zeros(1,order);
    for n = 1:order
        index_l_g = n^2;
        index_r_g = index_l_g+n;
        index_l_h = index_r_g+1;
        index_r_h = index_l_h+(n-1);
        gnm = GS(index_l_g:index_r_g)';
        hnm = [0,GS(index_l_h:index_r_h)'];
        Sn(n) = (n+1)*sum(gnm.^2+hnm.^2);
    end
end

function errorf(n,value1,value2,lable)
    subplot(2,2,n);
    plot(value1,'-o','LineWidth',1.5);
    hold on;
    plot(value2,'-p','LineWidth',1.5);
    xlabel('n of stations');
    ylabel('value(nT)');
    title(lable);
end