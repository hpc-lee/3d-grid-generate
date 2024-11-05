clear all;
close all;
clc;

% -------------------------- parameters input -------------------------- %
% file and path name
parfnm='../project/test.json';
output_dir='../project/output';

% which grid profile to plot
% subs=[30,1,1];    
% subc=[1,-1,-1];   % '-1' to plot all points in this dimension
% subt=[1,4,4];
subs=[1,150,1];    
subc=[-1,1,-1];   % '-1' to plot all points in this dimension
subt=[4,1,4];
% figure control parameters
flag_km     = 1;
flag_emlast = 1;
flag_print  = 1;
flag_title  = 0;
scl_daspect = [1 1 1];

%-----------------------------------------------------------
%-- load coord
%-----------------------------------------------------------
coordinfo=locate_coord(parfnm,output_dir,subs,subc,subt);
[x,y,z]=gather_coord(coordinfo,output_dir);
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
hid = figure;
set(hid,'BackingStore','on');
if subc(1) == 1
    plot(permute(y,[2,1]),permute(z,[2,1]),'k-');
    hold on;
    plot(y,z,'k-');
    xlabel(['y-axis (' str_unit ')']);
    ylabel(['z-axis (' str_unit ')']);
     
elseif subc(2) == 1
    plot(permute(x,[2,1]),permute(z,[2,1]),'k-');
    hold on;
    plot(x,z,'k-');
    xlabel(['x-axis (' str_unit ')']);
    ylabel(['z-axis (' str_unit ')']);
     
elseif subc(3) == 1
    plot(permute(x,[2,1]),permute(y,[2,1]),'k-');
    hold on;
    plot(x,y,'k-');
    xlabel(['x-axis (' str_unit ')']);
    ylabel(['y-axis (' str_unit ')']);
end
% h = surf(x,y,z);
set(gca,'layer','top');
set(gca,'FontSize',10,FontWeight='bold');
set(gcf,'color','white','renderer','painters');
set(gcf,'Position',[0,0,1920,1200]);
% set(h, 'FaceColor', 'none');
grid off;
box off;
% xlim([5.5,7.5]);
% ylim([0,1.2]);
axis equal tight;

% title
if flag_title
  if (subc(1) == 1)
    gridtitle='YOZ-Grid';
  elseif (subc(2) == 1)
    gridtitle='XOZ-Grid';
  elseif (subc(3) == 1)
    gridtitle='XOY-Grid';
  end
  title(gridtitle);
end

% save and print figure
if flag_print
  print(gcf,['grid2.png'],'-r300','-dpng');
end
