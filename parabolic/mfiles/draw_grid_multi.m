clear all;
close all;
clc;

% -------------------------- parameters input -------------------------- %
% file and path name
parfnm='../project1/test.json';
output_dir='../project1/output';

% profiles to plot
% profile 1
subs{1}=[401,1,1];      % start from index '1'
subc{1}=[1,-1,-1];     % '-1' to plot all points in this dimension
subt{1}=[1,1,1];
% profile 2
subs{2}=[1,1,1];      % start from index '1'
subc{2}=[-1,1,-1];     % '-1' to plot all points in this dimension
subt{2}=[1,1,1];
% profile 3
subs{3}=[1,1,201];      % start from index '1'
subc{3}=[-1,-1,1];     % '-1' to plot all points in this dimension
subt{3}=[1,1,1];
% profile 4
subs{4}=[201,1,1];      % start from index '1'
subc{4}=[1,-1,-1];     % '-1' to plot all points in this dimension
subt{4}=[1,1,1];
% profile 5
subs{5}=[1,149,1];      % start from index '1'
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
a=4;
h = figure;
set(h,'BackingStore','on');
surf(x{1}(149:a:301,1:a:201),y{1}(149:a:301,1:a:201),z{1}(149:a:301,1:a:201),'FaceColor', [1,1,1]);
hold on;
surf(x{2}(1:a:201,1:a:201),y{2}(1:a:201,1:a:201),z{2}(1:a:201,1:a:201),'FaceColor', [1,1,1]);
surf(x{3}(1:a:201,1:a:301),y{3}(1:a:201,1:a:301),z{3}(1:a:201,1:a:301),'FaceColor', [1,1,1]);
surf(x{3}(201:a:401,149:a:301),y{3}(201:a:401,149:a:301),z{3}(201:a:401,149:a:301),'FaceColor', [1,1,1]);
surf(x{4}(1:a:149,1:a:201),y{4}(1:a:149,1:a:201),z{4}(1:a:149,1:a:201),'FaceColor', [1,1,1]);
surf(x{5}(201:a:401,1:a:201),y{5}(201:a:401,1:a:201),z{5}(201:a:401,1:a:201),'FaceColor', [1,1,1]);
xlabel('X axis (km)',rotation=-20,position=[3.7 -3.7]);
ylabel('Y axis (km)',rotation=45,position=[6.0 -2.0]);
zlabel('Z axis (km)');
set(gca,'layer','top');
set(gca,'FontSize',15,FontWeight='bold');
set(gcf,'color','white');
set(gcf,'Position',[0,0,1920,1200]);
view(30,30)
axis equal tight;
box off;
grid off;
% cid=colorbar;
% cid.Label.String='km';
% set(cid,'Location','southoutside','Position',[0.21,0.038,0.61,0.025]);
% colormap(gray);

% save and print figure
if flag_print
  print(gcf,['grid1.png'],'-r400','-dpng');
end
