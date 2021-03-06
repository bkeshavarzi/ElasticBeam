clear
clc

fid=fopen('CSE_Node.txt','r');
line=fgetl(fid);
irow=0;

while (line >0)
    
    lin=strsplit(line,' ');
    irow=irow+1;
    NodeData(irow,1)=str2double(lin{2});
    NodeData(irow,2)=str2double(lin{3});
    NodeData(irow,3)=str2double(lin{4});
    line=fgetl(fid);
    
end

for i=1:size(NodeData,1)
    
   plot(NodeData(i,2),NodeData(i,3),'o')
   hold on
   
end

fclose(fid);

fid=fopen('CSE_Element.txt','r');
line=fgetl(fid);
irow=0;

while (line >0)
    
    lin=strsplit(line,' ');
    irow=irow+1;
    ElementData(irow,1)=str2double(lin{1});
    ElementData(irow,2)=str2double(lin{2});
    ElementData(irow,3)=str2double(lin{3});
    ElementData(irow,4)=str2double(lin{4});
    Cord_data(irow,1)=(1/3)*(NodeData(ElementData(irow,2),2)+NodeData(ElementData(irow,3),2)+NodeData(ElementData(irow,4),2));
    Cord_data(irow,2)=(1/3)*(NodeData(ElementData(irow,2),3)+NodeData(ElementData(irow,3),3)+NodeData(ElementData(irow,4),3));
    
    line=fgetl(fid);
    
end

for i=1:size(ElementData,1)
    
   plot([NodeData(ElementData(i,2),2),NodeData(ElementData(i,3),2)],[NodeData(ElementData(i,2),3),NodeData(ElementData(i,3),3)],'LineWidth',0.001)
   hold on
   plot([NodeData(ElementData(i,2),2),NodeData(ElementData(i,4),2)],[NodeData(ElementData(i,2),3),NodeData(ElementData(i,4),3)],'LineWidth',0.001)
   hold on
   plot([NodeData(ElementData(i,3),2),NodeData(ElementData(i,4),2)],[NodeData(ElementData(i,3),3),NodeData(ElementData(i,4),3)],'LineWidth',0.001)
   
end

fid=fopen('Displacement.txt');
line=fgetl(fid);
irow=0;

while (line >0)
    
lin=strsplit(line,'\t');
irow=irow+1;

DispData(irow,1)=str2num(lin{1})+1; %Node id
DispData(irow,2)=str2double(lin{2}); %x cord
DispData(irow,3)=str2double(lin{3}); %y cord
DispData(irow,4)=str2double(lin{4}); %u_x
DispData(irow,5)=str2double(lin{5}); %u_y

line=fgetl(fid);

end

close all

% figure()
% x_cord=(NodeData(:,2))';
% y_cord=(NodeData(:,3))';
% [xq,yq] = meshgrid(min(x_cord):0.1:max(x_cord), min(y_cord):0.1:max(y_cord));
% zq = griddata(x_cord,y_cord,DispData(:,4),xq,yq);
% contourf(xq,yq,zq);
% colorbar;



% figure()
% x_cord=(NodeData(:,2))';
% y_cord=(NodeData(:,3))';
% [xq,yq] = meshgrid(min(x_cord):0.1:max(x_cord), min(y_cord):0.1:max(y_cord));
% zq = griddata(x_cord,y_cord,DispData(:,5),xq,yq);
% contourf(xq,yq,zq);
% colorbar;
% grid on


fid=fopen('Stress_Strain.txt','r');
line=fgetl(fid);
irow=0;

while (line >0)
   
    irow=irow+1;
    lin=strsplit(line,'\t');
    
    SSData(irow,1)=str2double(lin{1}); %Epsilon xx
    SSData(irow,2)=str2double(lin{2}); %Epsilon yy
    SSData(irow,3)=str2double(lin{3}); %Epsilon xy
    
    SSData(irow,4)=str2double(lin{4}); %Sigma xx
    SSData(irow,5)=str2double(lin{5}); %Sigma yy
    SSData(irow,6)=str2double(lin{6}); %Sigma xy
    
    SSData(irow,7)=str2double(lin{7}); %PSigma xx
    SSData(irow,8)=str2double(lin{8}); %PSigma yy
    SSData(irow,9)=str2double(lin{9}); %PSigma xy
    
    SSData(irow,10)=str2double(lin{10}); %PStrain xx
    SSData(irow,11)=str2double(lin{11}); %PStrain yy
    SSData(irow,12)=str2double(lin{12}); %PStrain zz
    
    line=fgetl(fid);
    
end

for ielem=1:size(ElementData,1)
    
    inode=ElementData(ielem,2);
    jnode=ElementData(ielem,3);
    knode=ElementData(ielem,4);
    
    x_cord=[NodeData(inode,2),NodeData(jnode,2),NodeData(knode,2)];
    y_cord=[NodeData(inode,3),NodeData(jnode,3),NodeData(knode,3)];
    x_ave=(1/3)*(x_cord);
    y_ave=(1/3)*(y_cord);
    
    %u_x=[DispData(inode,4),DispData(jnode,4),DispData(knode,4)];
    patch(x_ave,y_ave,SSData(ielem,6),'LineStyle','none')
    
    hold on
    
end

colorbar;

% figure()
% [xq,yq] = meshgrid(min(Cord_data(:,1)):0.1:max(Cord_data(:,1)), min(Cord_data(:,2)):0.1:max(Cord_data(:,2)));
% zq = griddata(Cord_data(:,1),Cord_data(:,2),SSData(:,6),xq,yq);
% contourf(xq,yq,zq);
% colorbar;
