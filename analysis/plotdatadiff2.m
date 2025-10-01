function plotdatadiff2(file, traj, i, dwell, till, optionsave)
%function plotdatadiff2(file, traj, i, dwell, fractal,optionsave)
% creates one .tiff file per trajectory with diffusion analysis data
% Pc: averaged (column 7)
% Marianne Renner - 07/12 - SPTrack_v4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                              
if nargin<6
    optionsave=0;
end

till=till/1000;

f=figure;
if size(traj.fill,2)>4
   filling=traj.fill;
   index=find(filling(:,6)>0 | filling(:,6)<0);
   meanfillingsyn=mean(filling(index,5));
   control=1;
else
    control=0;
end

% plot 1: 
subplot(2,3,1) % m x n plots, p active
axis([0 10 0 8]);
text(1,7,['File : ',file]); 
text(1,6,['Trajectory # ',num2str(traj.coord(1,1))]);
text(1,5,[num2str(size(traj.coord,1)),' points']);
text(1,4,'D:'); text(4,4, sprintf('%6.5f',traj.D));
if control==1
    text(1,3,'MeanPc:'); text(6,3, sprintf('%6.2f',meanfillingsyn));
end
if isfield(traj,'dwell')
    if isempty(traj.dwell)==0
        text(1,2,'Dwell:'); text(4,2, sprintf('%4.2f',dwell));
        text(1,1,'Time syn:'); text(6,1, sprintf('%4.2f',traj.dwell(1)));
    else
        text(1,2,'Extrasynaptic'); 
    end
else
    % no loc
    text(1,2,'Not localized'); 
end
    
%plot 2:trajectory
subplot(2,2,2) % m x n plots, p active               
title('Trajectory')
for h=1:size(traj.coord,1)-1
  plot(traj.coord(h:h+1,3),traj.coord(h:h+1,4),'-b');   
  hold on
end

%  plot 3: MSD
subplot(2,2,3) % m x n plots, p active
if isfield(traj,'fourier') & isempty(traj.fourier)==0
    plot(traj.fourier(:,1),traj.fourier(:,3)), grid on
    [mp, indexfou]=max(traj.fourier(:,3));
    periodo=traj.fourier(indexfou,2);
    xlabel('proportion')
    title(['Per: ',num2str(periodo)])
else
    xlabel('t lag')
    ylabel('MSD')
    plot(traj.msd(:,1),traj.msd(:,2),'-b');
    msd=traj.msd;
    
   %save([file,'-',num2str(traj.coord(1,1)),'.msd'],'msd','-ascii');
    
end
% plot 4: filling
subplot(2,2,4) % m x n plots, p active
xlabel('frames')
ylabel('filling')

if size(traj.fill,2)>4
   plot(traj.fill(:,1),traj.fill(:,5),'-k');
%   plot(traj.fill(:,1),traj.fill(:,7),'-k');
   filling=[];
   filling=traj.fill;
   newPC=[];
   for j=1:size(filling,1)
       index=find(traj.coord(:,2)==filling(j,1));
       %frame
       if isempty(index)==0
           frame=filling(j,1);
           newPC(j,1)=frame;
           newPC(j,2)=frame*till; % frame ms
           newPC(j,3)=traj.coord(index(1),3); % x
           newPC(j,4)=traj.coord(index(1),4); % y
           newPC(j,5)=filling(j,5);
       end
   end
   if isempty(newPC)==0
       save([file,'-',num2str(traj.coord(1,1)),'.txt'],'newPC','-ascii');
   end
end

if optionsave==1
    % save 
    presentimage=getframe(f);  % gets the figure for the movie
    [image,Map] = frame2im(presentimage);
    imwrite(image,[file,'-',num2str(i),'.tif'],'tif');
else
    pause
end
close(f)
clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%