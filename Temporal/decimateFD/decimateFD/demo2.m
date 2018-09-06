m = 2;
t = 0:.00025:1;
x = sin(2*pi*30*t) + sin(2*pi*60*t);
subplot 211
stem(0:120,x(1:121),'filled','markersize',3);
ylabel 'Original'; grid on
subplot 212
yp = decimateFD(x,m);
stem(0:60,yp(1:61),'filled','markersize',3,'color','r');
xlabel 'Sample number',ylabel 'Decimated'; grid on;
print('-dpng','-painters','-r600','decimateFD-2.png');