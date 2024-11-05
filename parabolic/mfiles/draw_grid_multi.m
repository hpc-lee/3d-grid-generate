clear all;
close all;
clc;

% -------------------------- parameters input -------------------------- %
% file and path name
parfnm='../project/test.json';
output_dir='../project/output';

% profiles to plot
% profile 1
subs{1}=[1,1,301];      % start from index '1'
subc{1}=[-1,-1,1];     % '-1' to plot all points in this dimension
subt{1}=[1,1,1];
% profile 2
subs{2}=[241,1,1];      % start from index '1'
subc{2}=[1,-1,-1];     % '-1' to plot all points in this dimension
subt{2}=[1,1,1];
% profile 3
subs{3}=[541,1,1];      % start from index '1'
subc{3}=[1,-1,-1];     % '-1' to plot all points in this dimension
subt{3}=[1,1,1];
% profile 4
subs{4}=[1,1,1];      % start from index '1'
subc{4}=[-1,1,-1];     % '-1' to plot all points in this dimension
subt{4}=[1,1,1];
% profile 5
subs{5}=[1,241,1];      % start from index '1'
subc{5}=[-1,1,-1];     % '-1' to plot all points in this dimension
subt{5}=[1,1,1];

% figure control parameters
flag_km     = 1;
flag_emlast = 1;
flag_print  = 1;
flag_title  = 0;
scl_daspect = [1 1 1];


%-----------------------------------------------------------
%-- set figure
%-----------------------------------------------------------
hid = figure;
set(hid,'BackingStore','on');
% load data and plot
for i=1:length(subs)
  coordinfo{i}=locate_coord(parfnm,output_dir,subs{i},subc{i},subt{i});
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
end
a=6;
surf(x{1}(1:a:241,1:a:541),y{1}(1:a:241,1:a:541),z{1}(1:a:241,1:a:541));
hold on;
surf(x{1}(241:a:541,241:a:541),y{1}(241:a:541,241:a:541),z{1}(241:a:541,241:a:541));
surf(x{2}(1:a:241,1:a:301),y{2}(1:a:241,1:a:301),z{2}(1:a:241,1:a:301));
surf(x{3}(241:a:541,1:a:301),y{3}(241:a:541,1:a:301),z{3}(241:a:541,1:a:301));
surf(x{4}(1:a:241,1:a:301),y{4}(1:a:241,1:a:301),z{4}(1:a:241,1:a:301));
surf(x{5}(241:a:541,1:a:301),y{5}(241:a:541,1:a:301),z{5}(241:a:541,1:a:301));
%set(h, 'FaceColor', 'none');
xlabel('x-axis (km)',FontSize=15,FontWeight='bold',rotation=-50,position=[8 1.5]);
ylabel('y-axis (km)',FontSize=15,FontWeight='bold',rotation=10,position=[10.5 4.5]);
zlabel('z-axis (km)');
%xlabel('x-axis (km)',FontSize=15,FontWeight='bold',rotation=-10,position=[9 -5.0]);
%ylabel('y-axis (km)',FontSize=15,FontWeight='bold',rotation=20,position=[12.0 -2]);
%zlabel('z-axis (km)',FontSize=15,FontWeight='bold',rotation=0,position=[-3.5 9.5]);
set(gca,'layer','top');
set(gca,'FontSize',10,FontWeight='bold');
set(gcf,'color','white');
set(gcf,'Position',[0,0,1920,1200]);
%view(35,16)
view(60,30)
% xlim([5.5,7.5]);
% ylim([0,1.2]);
axis equal tight;
box off;
grid off;
cid=colorbar;
cid.Label.String='km';
set(cid,'Location','southoutside','Position',[0.21,0.038,0.61,0.025]);
% colormap(gray);

% save and print figure
if flag_print
  print(gcf,['grid1.png'],'-r300','-dpng');
end
