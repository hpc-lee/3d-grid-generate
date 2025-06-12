clear all;
close all;
clc;
addmypath
% -------------------------- parameters input -------------------------- %
% file and path name
parfnm='../project1/test.json';
output_dir='../project1/output';

% media profiles to plot
% profile 1
subs{1}=[201,1,1];      % start from index '1'
subc{1}=[1,-1,-1];     % '-1' to plot all points in this dimension
subt{1}=[1,4,4];
% profile 2
subs{2}=[1,201,1];      % start from index '1'
subc{2}=[-1,1,-1];     % '-1' to plot all points in this dimension
subt{2}=[4,1,4];

% variable to plot
% 'Vp', 'Vs', 'rho', 'lambda', 'mu'
varnm='Vs';

% figure control parameters
flag_km     = 1;
flag_print  = 1;
flag_clb    = 1;
flag_title  = 1;
clrmp       = 'parula';
% ---------------------------------------------------------------------- %

% figure plot
hid=figure;
set(hid,'BackingStore','on');

% load data and plot
for i=1:length(subs)
    
    % locate media
    mediainfo{i}=locate_media(parfnm,output_dir,subs{i},subc{i},subt{i});
    
    % get coordinate data
    [x{i},y{i},z{i}]=gather_coord(mediainfo{i},output_dir);
    %- set coord unit
    if flag_km
       x{i}=x{i}/1e3;
       y{i}=y{i}/1e3;
       z{i}=z{i}/1e3;
       str_unit='km';
    else
       str_unit='m';
    end
    
    % gather media
    switch varnm
        case 'Vp'
            rho=gather_media(mediainfo{i},'rho',output_dir);
               mu=gather_media(mediainfo{i},'mu',output_dir);
               lambda=gather_media(mediainfo{i},'lambda',output_dir);
               v{i}=( (lambda+2*mu)./rho ).^0.5;
            v{i}=v{i}/1e3;
        case 'Vs'
            rho=gather_media(mediainfo{i},'rho',output_dir);
            mu=gather_media(mediainfo{i},'mu',output_dir);
            v{i}=( mu./rho ).^0.5;
            v{i}=v{i}/1e3;
        case 'rho'
            v{i}=gather_media(mediainfo{i},varnm,output_dir);
            v{i}=v{i}/1e3;
    end
    
    % media show
    surf(x{i},y{i},z{i},v{i});
    hold on;
end
plot3(x{1}(:,401),y{1}(:,401),z{1}(:,401),'r',LineWidth=1.5);
hold on;
plot3(x{1}(1,:),y{1}(1,:),z{1}(1,:),'r',LineWidth=1.5);
plot3(x{2}(:,401),y{2}(:,401),z{2}(:,401),'r',LineWidth=1.5);
plot3(x{3}(1,:),y{3}(1,:),z{3}(1,:),'r',LineWidth=1.5);
plot3(x{3}(:,641),y{3}(:,641),z{3}(:,641),'r',LineWidth=1.5);
xlabel(['X axis (' str_unit ')'],FontSize=15,FontWeight='bold',rotation=-50,position=[12.5 0.0]);
ylabel(['Y axis (' str_unit ')'],FontSize=15,FontWeight='bold',rotation=5,position=[15.5 3.5]);
zlabel(['Z axis (' str_unit ')'],FontSize=15,FontWeight='bold',rotation=0,position=[-5.0 5.5]);
view(70,20);

set(gca,'layer','top');
%set(gcf,'color','white','renderer','painters');
set(gcf,'color','white');

% shading
shading interp;
%shading flat;
% colorbar range/scale
if exist('scl_caxis','var')
    caxis(scl_caxis);
end
axis equal tight;
% colormap and colorbar
if exist('clrmp')
    colormap(clrmp);
end
if flag_clb
    cid=colorbar;
    if strcmp(varnm,'Vp') || strcmp(varnm,'Vs')
        cid.Label.String='(km/s)';
    end
    if strcmp(varnm,'rho')
        cid.Label.String='g/cm^3';
    end
%     set(cid,'Location','southoutside','Position',[0.24,0.038,0.53,0.02]);
end
grid off
box on
set(gca, FontWeight='bold', FontSize=15);
set(gcf,'Position',[0,0,1920,1200]);
% title
if flag_title
%     title(varnm,FontSize=15,FontWeight='bold');
end

% save and print figure

if flag_print
   print(gcf,'3d-media.png','-r400','-dpng');
end


