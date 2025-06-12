clear all;
close all;
clc;

% -------------------------- parameters input -------------------------- %
% file and path name
parfnm='../project2/test.json';
output_dir='../project2/output';

% media profiles to plot
% profile 1
subs{1}=[101,1,1];      % start from index '1'
subc{1}=[1,-1,-1];     % '-1' to plot all points in this dimension
subt{1}=[1,10,10];
% profile 2
subs{2}=[1,1,1];      % start from index '1'
subc{2}=[-1,1,-1];     % '-1' to plot all points in this dimension
subt{2}=[10,1,10];
% profile 3
subs{3}=[1,1,51];      % start from index '1'
subc{3}=[-1,-1,1];     % '-1' to plot all points in this dimension
subt{3}=[10,10,1];

% figure control parameters
flag_km     = 1;
flag_print  = 1;
flag_clb    = 0;
flag_title  = 1;
clrmp       = 'parula';
% ---------------------------------------------------------------------- %

% figure plot
hid=figure;
set(hid,'BackingStore','on');

% load data and plot
for i=1:length(subs)
    
    % locate media
    coordinfo{i}=locate_coord(parfnm,output_dir,subs{i},subc{i},subt{i});
    
    % get coordinate data
    [x{i},y{i},z{i}]=gather_coord(coordinfo{i},output_dir);
    %- set coord unit
    if flag_km
       x{i}=x{i}/1e3;
       y{i}=y{i}/1e3;
       z{i}=z{i}/1e3;
       str_unit='km';
    else
       str_unit='m';
    end

    plot3(x{i},y{i},z{i},'k-',LineWidth=2);
    hold on;
    plot3(x{i}',y{i}',z{i}','k-',LineWidth=2);
end
xticks([]); 
yticks([]); 
zticks([]); 
% xlabel(['X axis (' str_unit ')'],FontSize=15,FontWeight='bold',rotation=-50,position=[12.5 0.0]);
% ylabel(['Y axis (' str_unit ')'],FontSize=15,FontWeight='bold',rotation=5,position=[15.5 3.5]);
% zlabel(['Z axis (' str_unit ')'],FontSize=15,FontWeight='bold',rotation=0,position=[-5.0 5.5]);
view(40,40);

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
%     cid=colorbar;
%     set(cid,'Location','southoutside','Position',[0.24,0.038,0.53,0.02]);
end
grid off
box off
axis off;
set(gca, FontWeight='bold', FontSize=15);
set(gcf,'Position',[0,0,1920,1200]);
% title
if flag_title
%     title(varnm,FontSize=15,FontWeight='bold');
end

% save and print figure

if flag_print
   print(gcf,'3d-grid2.png','-r400','-dpng');
end


