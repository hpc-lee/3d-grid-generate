clear all;
close all;
clc;
addmypath

% -------------------------- parameters input -------------------------- %
% file and path name
parfnm='../project/test.json';
output_dir='../project/output';

% % which grid profile to plot
% subs=[200,1,1];    
% subc=[1,-1,-1];   % '-1' to plot all points in this dimension
% subt=[1,1,1];
% subs=[1,1,100];    
% subc=[-1,-1,1];   % '-1' to plot all points in this dimension
% subt=[1,1,1];
subs=[1,150,1];    
subc=[-1,1,-1];   % '-1' to plot all points in this dimension
subt=[3,2,3];

% figure control parameters
flag_km     = 0;
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
varnm = 'orth_xizt';
% varnm = 'jacobi';
% varnm = 'step_zt';
% varnm = 'smooth_et';
%-----------------------------------------------------------
%-- load coord
%-----------------------------------------------------------

qualityinfo=locate_quality(parfnm,output_dir,subs,subc,subt);
[x,y,z]=gather_coord(qualityinfo,output_dir);

v=gather_quality(qualityinfo,output_dir,varnm);

%- set coord unit
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
% shading flat;
% colorbar range/scale
if exist('scl_caxis','var')
    caxis(scl_caxis);
end
% axis daspect
if exist('scl_daspect')
    daspect(scl_daspect);
end
axis equal
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
