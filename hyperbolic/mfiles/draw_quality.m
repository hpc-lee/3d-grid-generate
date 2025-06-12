clear all;
close all;
clc;

% -------------------------- parameters input -------------------------- %
% file and path name
parfnm='../project/test.json';
output_dir='../project/output';

% which grid profile to plot
%subs=[201,1,1];    
%subc=[1,-1,-1];   % '-1' to plot all points in this dimension
%subt=[1,1,1];
subs=[1,201,1];    
subc=[-1,1,-1];   % '-1' to plot all points in this dimension
subt=[1,1,1];

% figure control parameters
flag_km     = 1;
flag_emlast = 1;
flag_print  = 1;
flag_clb    = 1;

flag_title  = 1;
% scl_daspect = [1 1 1];
% scl_caxis = [1,1.1];
% scl_caxis = [40,90];
clrmp       = 'parula';

% varable to plot
% 'orth_xiet', 'orth_xizt',
% 'orth_etzt', 'jacobi',  
% 'smooth_xi', 'smooth_et','smooth_zt',
% % 'step_xi', 'step_et', 'step_zt'
% varnm = 'orth_xizt';
varnm='smooth_zt';
% varnm='smooth_et';
% varnm='smooth_xi';
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
    xlabel(['Y axis (' str_unit ')']);
    ylabel(['Z axis (' str_unit ')']);
     
elseif subc(2) == 1
    pcolor(x,z,v);
    xlabel(['X axis (' str_unit ')']);
    ylabel(['Z axis (' str_unit ')']);
     
elseif subc(3) == 1
    pcolor(x,y,v);
    xlabel(['X axis (' str_unit ')']);
    ylabel(['Z axis (' str_unit ')']);
end

shading interp;
set(gca,'layer','top');
% set(gcf,'color','white','renderer','painters');
set(gcf,'color','white');
set(gcf,'Position',[200,200,1200,800]);
set(gca,'FontSize',15,FontWeight='bold');

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
  set(get(cid,'Title'),'string','degree');
end
if flag_title
%       title('J',FontSize=15);
%     title('Q\_\xi\eta',FontSize=15);
    title('Q\_\xi\zeta',FontSize=15);
%     title('Q\_\eta\zeta',FontSize=15);
%     title('S\_\xi',FontSize=15);
%     title('S\_\eta',FontSize=15);
%     title('S\_\zeta',FontSize=15);
end

% save and print figure
if flag_print
  print(gcf,['03-3d-model2-',varnm,'.png'],'-r400','-dpng');
end
