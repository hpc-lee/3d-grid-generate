clear all;
close all;
clc;

% -------------------------- parameters input -------------------------- %
% file and path name
parfnm='../project1/test.json';
output_dir='../project1/output';

% which grid profile to plot
subs=[2,1,1];    
subc=[1,-1,-1];   % '-1' to plot all points in this dimension
subt=[1,1,1];

% figure control parameters
flag_km     = 1;
flag_emlast = 1;
flag_print  = 1;
flag_clb    = 1;

flag_title  = 1;
scl_daspect = [1 1 1];
% scl_caxis = [1,1.1];
clrmp       = 'parula';

% varable to plot
% 'orth_xiet', 'orth_xizt',
% 'orth_etzt', 'jacobi',  
% 'smooth_xi', 'smooth_et','smooth_zt',
% 'step_xi', 'step_et', 'step_zt'
varnm = 'orth_etzt';
% varnm='smooth_zt';
% varnm = 'jacobi';
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


shading interp;
set(gca,'layer','top');
set(gca,'FontSize',10,FontWeight='bold');
% set(gcf,'color','white','renderer','painters');
set(gcf,'color','white');
set(gcf,'Position',[200,200,650,400]);

% colorbar range/scale
if exist('scl_caxis','var')
    caxis(scl_caxis);
end
% axis daspect
if exist('scl_daspect')
    daspect(scl_daspect);
end
axis equal tight;
% colormap and colorbar
if exist('clrmp')
    colormap(clrmp);
end
if flag_clb
    cid=colorbar;
%   set(get(cid,'Title'),'string','degree');
end
if flag_title
%     title('Orthogonality\_\eta\zeta',FontSize=15);
%     title('Smooth\_\xi',FontSize=15);
%     title('Smooth\_\eta',FontSize=15);
    title('Smooth\_\zeta',FontSize=15);
end

% save and print figure
if flag_print
  print(gcf,[varnm,'.png'],'-r400','-dpng');
end
