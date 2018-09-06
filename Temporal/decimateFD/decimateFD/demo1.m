clear; clc; close all;

% load original signal sampled at 240 Hz
load original.mat
dt = 1 / 240;

n = 4; % downsampling factor

timespan = 0:dt:(length(original)-1)*dt; 

figure(1)
fsz = 16;
set(gcf,'position',[300 283 600 600]);
set(gca,'TickLength',[.0025 .0025]);

subplot(211)
plot(timespan,original);
title('Original signal','fontsize',fsz,'fontname','times');
set(gca,'fontname','times','fontsize',fsz);

subplot(212)
yp = downsample(original,n); 
timespan = 0:dt*n:(length(original)-1)*dt; 
plot(timespan,yp,'r'); hold on; 

yp = decimateFD(original,n); 
plot(timespan,yp,'g'); 
title('Decimated signal by 4 (MatLAB downsample vs decimateFD)','fontsize',fsz,'fontname','times');
legend('MatLAB downsample','decimateFD'); 
xlabel('Time (s)','FontSize',[fsz],'fontname','times');
set(gca,'fontname','times','fontsize',fsz); 

print('-dpng','-painters','-r600','decimateFD-1.png')