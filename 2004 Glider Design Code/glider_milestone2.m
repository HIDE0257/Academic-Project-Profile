%% weight at take off 
w = .295*9.81;

A_fuse = .139;
A_wing = .177;
A_can = .03;

w_fuse = (A_fuse)*w;
m_fuse = w_fuse/9.81;

w_wing = (A_wing)*w;
m_wing = w_wing/9.81;

w_can = (A_can)*w;
m_can = w_can/9.81; 

w_payload = .16;
m_payload = w_payload/9.81; 

w_to = w_fuse+w_wing+w_can+w_payload;
m_to = w_to/9.81;
%% center of gravity 
x_fuse = (.343/2)+((.186+.037+.062)/4);
x_wing = .5*(.103+.077);
x_can = .35/2;
x_payload = (1/12)*m_payload*(.039^2+.0698^2);

cg = (((x_wing*m_wing)+(x_can*m_can)+(x_payload*m_payload)+(x_fuse*m_fuse))/m_to);

mcg_fuse = w_fuse*x_fuse;
mcg_wing = w_wing*x_wing;
mcg_can = w_can*x_can;

mcg_payload = w_payload*x_payload;

Mcg = mcg_fuse+mcg_wing+mcg_can+mcg_payload;

%% wing loading 
S = A_wing + A_can+A_fuse;
wing_loading = w_to/S;

%% optimize horizontl tail 

c = .18;
ch = linspace(.01,1,100);
bh = linspace(.01,1,100);
sh = ch.*bh;

cv = linspace(0.01,1,100);
bv = linspace(0.01,1,100);
sv = cv.*bv;

vol_h = (((ch.*bh).*abs(((.25*ch)-cg)))./(0.1*c));
 

figure(1)
hold on
yyaxis left
plot(vol_h, ch);
title('Surface Area of Horizontal Tail to Acieve Reasonable V_H');
xline([0.3,0.6],'--',{'Min','Max'});
xlabel('Volume Coefficient for Horizontal Tail');
ylabel('Chord length (m)');
axis([0 0.8 0 0.7]);

yyaxis right 
plot(vol_h, sh);
ylabel('Surface Area (m^2)');
grid on
hold off

% figure(2)
% hold on
% yyaxis left
% plot(vol_v, cv);
% title('Surface Area of Horizontal Tail to Acieve Reasonable V_H');
% xline([0.3,0.6],'--',{'Min','Max'});
% xlabel('Volume Coefficient for Horizontal Tail');
% ylabel('Chord length (m)');
% axis([0 0.8 0 0.7]);
% 
% yyaxis right 
% plot(vol_v, sv);
% ylabel('Surface Area (m^2)');
% grid on
% hold off
