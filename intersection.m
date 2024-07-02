% clc
% array=[515 525 561 600 632 700 761 800 900 1000 1014 1750;
%     0 150 300 394 450 540 600 631 696 745 750 865];
% x=linspace(array(1,1),array(1,end),101)
% y=interp1(array(1,:),array(2,:),x,'pchip')
% x=transpose(x)
% y=transpose(y)
% %
% f=fit(y,x,'pchip')
% a=coeffvalues(f)
% plot(f,y,x)
% hold on
% % Equation of line that pass through origin
% x1=0:1000;
% slope=tan(51.5*pi/180);
% y1=slope*x1
% plot(x1,y1)

%% 

array=[515 525 561 600 632 700 761 800 900 1000 1014 1750;
         0 150 300 394 450 540 600 631 696  745  750  865];
x = array([2 1],:)'
f = fit(x(:,1),x(:,2),'cubicinterp')
df = fit(x(:,1),differentiate(f,x(:,1)),'cubicinterp');
xx = fzero(@(x)f(x)/x - df(x),[1 750]);
x1 = linspace(x(1,1),x(end,1),300);
plot(x1,f(x1),x1,x1*df(xx),xx,f(xx),'ro'); 

%% 

% x = -4 : .1 : 4;
% y1 = x + 3;
% y2 = x.^2 -4; 
% 
% plot(x,y1,'b', x,y2,'r')
% grid on
% title('Intersection')
% xlabel('x')
% ylabel('y') 
% 
% d = y1 ./ (y2 + eps);
% ix = find(d > .95 & d < 1.05);
% x_sol = x(ix)
% y1_sol = y1(ix)
% y2_sol = y2(ix) 
