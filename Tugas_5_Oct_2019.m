%% TUGAS FLUID SUBTITION                    NISFU ANUGERAH /101116027/ TEKNIK GEOFISIKA
% (A) --> simbol Interval A
% (B) --> simbol Interval B
% (1) --> Kondisi Awal
% (2) --> Setalah Terssubtitusi Fluida Lain

clear, clc, close all

data = load('Data_Log.txt');
depth = data(:,1);
rho = data(:,2);
Vp = data(:,3);
Vs = data(:,4);
So = data(:,5);
%KONDISI AWAL
SoA = 0.7; 
SoB = 0.6;
SwA = 1-SoA;
SwB = 1-SoB;
%POROSITAS
PorA = 0.3;
PorB = 0.35;

% Kasus
Swa = [(ones(size(SoA)))]; % untuk interval A
Swb = [(ones(size(SoB)))]; % untuk interval B

% Fraction
%% Interval A
FqA = (1-PorA)*0.7;
FcA = (1-PorA)*0.3;
FoA = 0.3*SoA;
FwA = 0.3*SwA;
%% Interval B 
FqB = (1-PorB)*0.6;
FcB = (1-PorB)*0.4;
FoB = PorB*SoB;
FwB = PorB*SwB;

% Bulk Moduli for each minerals and fluid
Kq = 36;
Kc = 75;
Ko = 1.6;
Kw = 2.38;
%% rho for each minerals and fluids
rhoq = 2.65;
rhoc = 2.71;
rhoil=0.8;
rhow = 1;

%%Voight (upper-bound)
KupA = FqA*Kq + FcA*Kc + FoA*Ko + FwA*Kw; 
KupB = FqB*Kq + FcB*Kc + FoB*Ko + FwB*Kw;
%%Reuss (lower-bound)
KlowA = 1/((FqA/Kq + FcA/Kc + FoA/Ko + FwA/Kw));
KlowB = 1/((FqB/Kq + FcB/Kc + FoB/Ko + FwB/Kw));
%%Voight-Reuss-Hill for Kmineral
KmA = (KupA + KlowA')/2;
KmB = (KupB + KlowB')/2;

%% Interval A
depthA = depth(1:34);
rhoA1 = rho(1:34);%interval
VpA1 = Vp(1:34); %Vinterval (m/s)
VsA1 = Vs(1:34); %Vinterval (m/s)
KsatA1 = rhoA1.*((VpA1.*10^-3).^2 - (4/3*(VsA1.*10^-3).^2)); % Using Interval (*10^-3 agar km/s)
MsatA1 = rhoA1.*((VsA1.*10^-3).^2);
KflA1 = 1./(SoA./Ko + SwA./Kw);
KflA2 = 1./(Swa./Kw);
A = ((KsatA1./(KmA-KsatA1)- KflA1./(PorA.*(KmA-KflA1)) + KflA2./(PorA.*(KmA-KflA2))));
KsatA2= (A.*KmA)./(A+1);
rhoA2 = (rhoA1 + PorA*(0.2)).*ones(size(rhoA1));
MsatA2 = MsatA1;

VpA2 = (sqrt((KsatA2+(4/3*MsatA2))./rhoA2)).*10^3; %m/s
VsA2 = (sqrt(MsatA2./rhoA2)).*10^3; %m/s

%% Interval B
depthB = depth(35:60);
rhoB1 = rho(35:60);%interval
VpB1 = Vp(35:60); %Vinterval (m/s)
VsB1 = Vs(35:60); %Vinterval (m/s)
KsatB1 = rhoB1.*((VpB1.*10^-3).^2 - (4/3*(VsB1.*10^-3).^2)); % Using Interval (*10^-3 agar km/s)
MsatB1 = rhoB1.*((VsB1.*10^-3).^2);
%%Fluid Subtitution
KflB1 = 1./(SoB./Ko + SwB./Kw);
KflB2 = 1./(Swb./Kw);
B = ((KsatB1./(KmB-KsatB1)- KflB1./(PorB.*(KmB-KflB1)) + KflB2./(PorB.*(KmB-KflB2))));
KsatB2= (B.*KmA)./(B+1);
rhoB2 = (rhoB1 + PorB*(0.2)).*ones(size(rhoB1));
MsatB2 = MsatB1;

VpB2 = (sqrt((KsatB2+(4/3*MsatB2))./rhoB2)).*10^3; %m/s
VsB2 = (sqrt(MsatB2./rhoB2)).*10^3; %m/s

%% selisih dari sesudah - sebelum
VPresA = (VpA2-VpA1);
VPresB = (VpB2-VpB1);

RHOresA = (rhoA2-rhoA1);
RHOresB = (rhoB2-rhoB1);

KsatresA =(KsatA2-KsatA1);
KsatresB =(KsatB2-KsatB1);


figure(1)
subplot(2,1,1);
plot(VpA1,depthA), hold on, plot(VpA2,depthA,'r')
xlabel({'VP','(m/s)'})
ylabel({'Depth','(m)'})
legend({'VpA1','VpA2'},'location','southeast')
legend('boxoff')
title('Interval A')
set(gca,'Ydir','reverse')
subplot(2,1,2);
plot(VpB1,depthB), hold on, plot(VpB2,depthB,'r')
xlabel({'VP','(m/s)'})
ylabel({'Depth','(m)'})
legend({'VpB1','VpB2'},'location','southeast')
legend('boxoff')
title('Interval B')
set(gca,'Ydir','reverse')

figure(2)
subplot(2,1,1);
plot(rhoA1,depthA), hold on, plot(rhoA2,depthA,'r')
xlabel({'Density','(g/cm3)'})
ylabel({'Depth','(m)'})
legend({'rhoA1','rhoA2'},'location','northeast')
legend('boxoff')
title('Interval A')
set(gca,'Ydir','reverse')
subplot(2,1,2);
plot(rhoB1,depthB), hold on, plot(rhoB2,depthB,'r')
xlabel({'Density','(g/cm3)'})
ylabel({'Depth','(m)'})
legend({'rhoB1','rhoB2'},'location','southeast')
legend('boxoff')
title('Interval B')
set(gca,'Ydir','reverse')

figure(3)
subplot(2,1,1);
plot(VPresA,depthA,'k')
xlabel({'VP','(m/s)'})
ylabel({'Depth','(m)'})
legend({'VpresA'},'location','southeast')
legend('boxoff')
title('Interval A')
set(gca,'Ydir','reverse')
subplot(2,1,2);
plot(VPresB,depthB,'k')
xlabel({'VP','(m/s)'})
ylabel({'Depth','(m)'})
legend({'VPresB'},'location','southeast')
legend('boxoff')
title('Interval B')
set(gca,'Ydir','reverse')

figure(4)
subplot(2,1,1);
plot(RHOresA,depthA,'r')
xlabel({'Density','(g/cm3)'})
ylabel({'Depth','(m)'})
legend({'RHOresA'},'location','southeast')
legend('boxoff')
title('Interval A')
set(gca,'Ydir','reverse')
subplot(2,1,2);
plot(RHOresB,depthB,'r')
xlabel({'Density','(g/cm3)'})
ylabel({'Depth','(m)'})
legend({'RHOresB'},'location','southeast')
legend('boxoff')
title('Interval B')
set(gca,'Ydir','reverse')

figure(5)
subplot(2,1,1);
plot(KsatA1,depthA), hold on, plot(KsatA2,depthA,'r')
xlabel({'Ksat','(GPa)'})
ylabel({'Depth','(m)'})
legend({'KsatA1','KsatA2'},'location','southeast')
legend('boxoff')
title('Interval A')
set(gca,'Ydir','reverse')
subplot(2,1,2);
plot(KsatB1,depthB), hold on, plot(KsatB2,depthB,'r')
xlabel({'Ksat','(GPa)'})
ylabel({'Depth','(m)'})
legend({'KsatB1','KsatB2'},'location','southeast')
legend('boxoff')
title('Interval B')
set(gca,'Ydir','reverse')

figure(6)
subplot(2,1,1);
plot(KsatresA,depthA)
xlabel({'Ksat','(GPa)'})
ylabel({'Depth','(m)'})
legend({'KsatresA'},'location','southeast')
legend('boxoff')
title('Interval A')
set(gca,'Ydir','reverse')
subplot(2,1,2);
plot(KsatresB,depthB)
xlabel({'Ksat','(GPa)'})
ylabel({'Depth','(m)'})
legend({'KsatresB'},'location','southeast')
legend('boxoff')
title('Interval B')
set(gca,'Ydir','reverse')

 




















