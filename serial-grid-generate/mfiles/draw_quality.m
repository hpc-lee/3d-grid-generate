clear all;
close all;
clc;
addmypath

% -------------------------- parameters input -------------------------- %
% file and path name
output_dir='../project/output';

% get nx ny and nz
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
subs=[200,1,1];     % index 1:nx 1:ny 1:nz
subc=[1,-1,-1];   % '-1' to plot all points in this dimension
subt=[1,1,1];
% subs=[1,100,1];     % index 1:nx 1:ny 1:nz
% subc=[-1,1,-1];   % '-1' to plot all points in this dimension
% subt=[1,1,1];
% subs=[1,1,100];     % index 1:nx 1:ny 1:nz
% subc=[-1,-1,1];   % '-1' to plot all points in this dimension
% subt=[1,1,1];

% figure control parameters
flag_km     = 1;
flag_emlast = 1;
flag_print  = 0;
flag_clb    = 1;
flag_title  = 1;
scl_daspect = [1 1 1];
clrmp       = 'parula';

% varable to plot 
% 'orth_xiet', 'orth_xizt',
% 'orth_etzt', 'jacobi',  
% 'smooth_xi', 'smooth_et','smooth_zt',
% 'step_xi', 'step_et', 'step_zt'
% varnm = 'step_zt';
varnm = 'orth_etzt';
%-----------------------------------------------------------
%-- load coord
%-----------------------------------------------------------

[x,y,z]=gather_coord(output_dir,subs,subc,subt);

v=gather_quality(output_dir,varnm,subs,subc,subt);

x = squeeze(x);
y = squeeze(y);
z = squeeze(z);
v = squeeze(v);

%- set coord unit
flag_km = 0;
if flag_km
   x=x/1e3;
   y=y/1e3;
   z=z/1e3;
   str_unit='km';
else
   str_unit='m';
end

%-----------------------------------------------------------
%-- set figure
%-----------------------------------------------------------
% figure plot
hid=figure;
set(hid,'BackingStore','on');

if subc(1) == 1
   pcolor(y,z,v);
   xlabel(['Y axis (' str_unit ')']);
   ylabel(['Z axis (' str_unit ')']);
   
elseif subc(2) == 1
   pcolor(x,z,v);
   xlabel(['X axis (' str_unit ')']);
   ylabel(['Z axis (' str_unit ')']);
   
elseif subc(3) == 1
   pcolor(x,y,v);
   xlabel(['X axis (' str_unit ')']);
   ylabel(['Y axis (' str_unit ')']);
end

set(gca,'layer','top');
set(gcf,'color','white','renderer','painters');

% shading
% shading interp;
%shading flat;
% colorbar range/scale
if exist('scl_caxis')
    caxis(scl_caxis);
end
% axis daspect
if exist('scl_daspect')
    daspect(scl_daspect);
end
axis tight
% colormap and colorbar
if exist('clrmp')
    colormap(clrmp);
end
if flag_clb
    cid=colorbar;
end

% title
if flag_title
    title(varnm,'interpreter','none');
end

% save and print figure
if flag_print
    width= 500;
    height=500;
    set(gcf,'paperpositionmode','manual');
    set(gcf,'paperunits','points');
    set(gcf,'papersize',[width,height]);
    set(gcf,'paperposition',[0,0,width,height]);
    print(gcf,[varnm '.png'],'-dpng');
end