clear all;
close all;
clc;

% -------------------------- parameters input -------------------------- %
% file and path name
parfnm='../project1/test.json';
output_dir='../project1/output';

% which grid profile to plot
% subs=[1,1,1];    
% subc=[-1,-1,1];   % '-1' to plot all points in this dimension
% subt=[2,2,2];

subs=[1,200,1];    
subc=[-1,1,-1];   % '-1' to plot all points in this dimension
subt=[4,1,4];

% figure control parameters
flag_km     = 0;
flag_print  = 0;
flag_clb    = 1;
flag_title  = 1;
clrmp       = 'parula';

% varable to plot 
% 'orth_xiet', 'orth_xizt',
% 'orth_etzt', 'jacobi',  
% 'smooth_xi', 'smooth_et','smooth_zt',
% 'step_xi', 'step_et', 'step_zt'
varnm = 'orth_xizt';
% varnm = 'step_xi';
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
    xlabel(['y-axis (' str_unit ')']);
    ylabel(['z-axis (' str_unit ')']);
     
elseif subc(2) == 1
    pcolor(x,z,v);
    xlabel(['x-axis (' str_unit ')']);
    ylabel(['z-axis (' str_unit ')']);
     
elseif subc(3) == 1
    pcolor(x,y,v);
    xlabel(['x-axis (' str_unit ')']);
    ylabel(['y-axis (' str_unit ')']);
end
     

set(gca,'layer','top');
set(gcf,'color','white');

% shading
shading interp;
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
end
% title
if flag_title
    title(varnm,'Interpreter','none');
end

% save and print figure
if flag_print
    print(gcf,[varnm '.png'],'-dpng');
end
