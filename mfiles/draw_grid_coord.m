clear all;
close all;
clc;
addmypath

% -------------------------- parameters input -------------------------- %
% file and path name
output_dir='../project/output';

% get nx and nz
fnm_coord=[output_dir,'/','coord.nc'];
if ~ exist(fnm_coord,'file')
   error([mfilename ': file ' fnm_coord 'does not exist']);
end

xyzc = nc_attget(fnm_coord,nc_global,'number_of_points');
xyzc = double(xyzc);
nx = xyzc(1);
ny = xyzc(2);
nz = xyzc(3);

% which grid profile to plot
subs=[1,100,1];     % index 1:nx 1:ny 1:nz
subc=[-1,1,-1];   % '-1' to plot all points in this dimension
subt=[1,1,1];

% figure control parameters
flag_km     = 1;
flag_emlast = 1;
flag_print  = 0;
flag_title  = 1;
scl_daspect = [1 1 1];
%-----------------------------------------------------------
%-- load coord
%-----------------------------------------------------------

[x,y,z]=gather_coord(output_dir,subs,subc,subt);

%- set coord unit
flag_km     = 0;
if flag_km
   x=x/1e3;
   y=y/1e3;
   z=z/1e3;
   str_unit='km';
else
   str_unit='m';
end

% plot
x=squeeze(x);
y=squeeze(y);
z=squeeze(z);
figure(1);
if subc(1) == 1
    plot(permute(y,[2,1]),permute(z,[2,1]),'k-');
    hold on;
    plot(y,z,'k-');
    xlabel(['Y axis (' str_unit ')']);
    ylabel(['Z axis (' str_unit ')']);
     
elseif subc(2) == 1
    plot(permute(x,[2,1]),permute(z,[2,1]),'k-');
    hold on;
    plot(x,z,'k-');
    xlabel(['X axis (' str_unit ')']);
    ylabel(['Z axis (' str_unit ')']);
     
elseif subc(3) == 1
    plot(permute(x,[2,1]),permute(y,[2,1]),'k-');
    hold on;
    plot(x,y,'k-');
    xlabel(['X axis (' str_unit ')']);
    ylabel(['Y axis (' str_unit ')']);
end

set(gca,'layer','top');
set(gcf,'color','white','renderer','painters');

% axis daspect
if exist('scl_daspect')
    daspect(scl_daspect);
end
axis tight;

% title
if flag_title
    if subc(1) == 1
        gridtitle='YOZ-Grid';
    elseif subc(2) == 1
        gridtitle='XOZ-Grid';
    elseif subc(3) == 1
        gridtitle='XOY-Grid';
    end
    title(gridtitle);
end

% save and print figure
if flag_print
    width= 500;
    height=500;
    set(gcf,'paperpositionmode','manual');
    set(gcf,'paperunits','points');
    set(gcf,'papersize',[width,height]);
    set(gcf,'paperposition',[0,0,width,height]);
    print(gcf,[gridtitle '.png'],'-dpng');
end
