clc
clear all

c = 30.4; %mm, chord
rho = 1.18; %kg/m^3, air density
mu = 1.846e-5;
TF = 1.2512;
patm = 100067.79;

figure
count2 = 1;
count3 = -6;
p2posavg = zeros(44,1);
alpha = [-6:24];
while count2<=31
    filechange = '2-5e5pressures/alpha_%d.dat';
    pRe2{count2}=load(sprintf(filechange,count3));
    count3 = count3+1;
    count = 1;
while count<=44
    pinf2pos = pRe2{count2}(:,4).*47.88;
    vinf2 = sqrt(2*pinf2pos/rho);
    p2pos = pRe2{count2}(count,5:14);
    p2pos = mean(p2pos)*47.88;
    p2posavg(count) = p2pos;
    count = count+1;
end
Cp2pos = (p2posavg)./(.5*rho*vinf2.^2);
X = pRe2{count2}(:,2);
Y = (Cp2pos);
intCpl2 = Y(1:23,1);
intCpu2 = Y(24:44,1);
Cl_2(count2) = trapz(X,Cp2pos)./c;
err25(count2) = std(Cl_2)./31;
plot(X./c,Y,'m')
set(gca,'Ydir','reverse');
title('Cp vs x/c, Re2.5e5')
xlabel('x/c')
ylabel('Cp')
grid on
hold on
count2 = count2+1;
end 

xf25 = importdata('red-airfoil_25e5.pol');
xfalpha25 = xf25(:,1);
xfCd25 = xf25(:,3);
xfCl25 = xf25(:,2);
figure
plot(alpha,Cl_2,'g--*',xfalpha25,xfCl25,'b--+')
title('Cl vs alpha, Re2.5e5')
xlabel('alpha')
ylabel('Cl')


figure
errorbar(alpha,Cl_2,err25,'r-')
title('Cl vs alpha, Re2.5e5')
xlabel('alpha')
ylabel('Cl')


figure
count2 = 1;
count3 = -6;
alphaD = [-6:3:24]';
Pdiffsave = zeros(200,11);
while count2<=11
    filechange = '25e5dragdata/alpha_%d.dat';
    pRe2drag{count2}=sortrows(load(sprintf(filechange,count3)));
    count3 = count3+3;
    count = 1;
while count<=200
    p2s = pRe2drag{count2}(count,3:12);
    p2savg = mean(p2s)*47.88;
    p2wake(count,1) = p2savg+patm;
    count = count+1;
end
    qinf2 = (pRe2drag{count2}(1,2).*47.88);
    vinf2 = sqrt(2*qinf2/rho);
    Pfstotal = qinf2+patm;
    Y = pRe2drag{1}(:,1);
    
    Pdiff = (Pfstotal-p2wake)./qinf2;
    Pdiffsave(:,count2) = Pdiff.*qinf2;
    Pdiffsave(Pdiffsave<0)=0;
    
    Cdcorr(count2) = -trapz(Pdiffsave(:,count2),Y).*1/(c*qinf2);
    errcd25(count2) = std(Cdcorr)./11;
    plot(Pdiff,Y./c)
    title('Pressure differential across wake, Re2.5e5')
    xlabel('Pressure differential')
    ylabel('Wake-rake location/c')
    hold on
count2 = count2+1;
end
hold off

figure
errorbar(alphaD,Cdcorr,errcd25,'r-')
title('Cd vs alpha, Re2.5e5')
xlabel('alpha')
ylabel('Cl')

figure
plot(alphaD,Cdcorr,'b--x',xfalpha25,xfCd25,'k--o')
title('Cd vs alpha')
xlabel('alpha')
ylabel('Cd')

figure
plot(Cdcorr,Cl_2(1:3:31),'b--x',xfCd25,xfCl25,'k--o')
title('Cd vs Cl, Re2.5e5')
xlabel('Cd')
ylabel('Cl')

%%
figure
count2 = 1;
count3 = -6;
p5posavg = zeros(44,1);
while count2<=31
    
    filechange = '5e5pressures/alpha_%d.dat';
    pRe5{count2}=load(sprintf(filechange,count3));
    
    count3 = count3+1;
    count = 1;
while count<=44
    
    pinf5pos = pRe5{count2}(:,4).*47.88;
    vinf5 = sqrt(2*pinf5pos/rho);
    p5pos = pRe5{count2}(count,5:14);
    p5pos = mean(p5pos)*47.88;
    p5posavg(count) = p5pos;

    count = count+1;
end

Cp5pos = (p5posavg)./(.5*rho*vinf5.^2);
X_5 = pRe5{count2}(:,2);
Y_5 = smooth(Cp5pos);
Cl_5(count2) = trapz(X_5,Cp5pos)./c;
err5(count2) = std(Cl_5)./31;
plot(X_5./c,Y_5,'b')
set(gca,'Ydir','reverse');
title('Cp vs x/c, Re5e5')
xlabel('x/c')
ylabel('Cp')
grid on
hold on

count2 = count2+1;
end 

xf5 = importdata('red-airfoil_5e5.pol');
xfalpha5 = xf5(:,1);
xfCd5 = xf5(:,3);
xfCl5 = xf5(:,2);
figure
plot(alpha,Cl_5,'g--*',xfalpha5,xfCl5,'b--+')
title('Cl vs alpha, Re5e5')
xlabel('alpha')
ylabel('Cl')


figure
errorbar(alpha,Cl_5,err5,'r-')
title('Cl vs alpha, Re5e5')
xlabel('alpha')
ylabel('Cl')

figure

count2 = 1;
count3 = -6;
alphaD = [-6:3:24]';
Pdiffsave = zeros(200,11);

while count2<=11
    filechange = '5e5dragdata/alpha_%d.dat';
    pRe5drag{count2}=sortrows(load(sprintf(filechange,count3)));
    count3 = count3+3;
    count = 1;
while count<=200
    p5s = pRe5drag{count2}(count,3:12);
    p5savg = mean(p5s)*47.88;
    p5wake(count,1) = p5savg+patm;
    count = count+1;
end
    qinf5 = (pRe5drag{count2}(1,2).*47.88);
    vinf5 = sqrt(2*qinf5/rho);
    Pfstotal = qinf5+patm;
    Y = pRe5drag{1}(:,1);
    
    Pdiff = (Pfstotal-p5wake)./qinf5;
    Pdiffsave(:,count2) = Pdiff.*qinf5;
    Pdiffsave(Pdiffsave<0)=0;
    
    Cdcorr5(count2) = -trapz(Pdiffsave(:,count2),Y).*1/(c*qinf5);
    errcd5(count2) = std(Cdcorr5)./11;
    plot(Pdiff,Y./c)
        title('Pressure differential across wake, Re5e5')
    xlabel('Pressure differential')
    ylabel('Wake-rake location/c')
    hold on
count2 = count2+1;
end
hold off

figure
errorbar(alphaD,Cdcorr5,errcd5,'r-')
title('Cd vs alpha, Re5e5')
xlabel('alpha')
ylabel('Cl')

figure
plot(alphaD,Cdcorr5,'b--x',xfalpha5,xfCd5,'k--o')
title('Cd vs alpha, Re5e5')
xlabel('alpha')
ylabel('Cd')

figure
plot(Cdcorr5,Cl_5(1:3:31),'b--x',xfCd5,xfCl5,'k--o')
title('Cd vs Cl, Re5e5')
xlabel('Cd')
ylabel('Cl')
%%
figure
count2 = 1;
count3 = -6;
p6posavg = zeros(44,1);
while count2<=31
    filechange = '6e5pressures/alpha_%d.dat';
    pRe6{count2}=load(sprintf(filechange,count3));
    count3 = count3+1;
    count = 1;
while count<=44
    pinf6pos = pRe6{count2}(:,4).*47.88;
    vinf6 = sqrt(2*pinf6pos/rho);
    p6pos = pRe6{count2}(count,5:14);
    p6pos = mean(p6pos)*47.88;
    p6posavg(count) = p6pos;
    count = count+1;
end
Cp6pos = (p6posavg)./(.5*rho*vinf6.^2);
X_6 = pRe6{count2}(:,2);
Y_6 = smooth(Cp6pos);
plot(X_6./c,Y_6,'g')
title('Cp vs x/c, Re6e5')
xlabel('x/c')
ylabel('Cp')
intCpl6 = Y_6(1:23,1);
intCpu6 = Y_6(24:44,1);
Cl_6(count2) = trapz(X_6,Cp6pos)./c;
err6(count2) = std(Cl_6)./31;
set(gca,'Ydir','reverse');
grid on
hold on
count2 = count2+1;
end 
xf6 = importdata('red-airfoil_6e5.pol');
xfalpha6 = xf25(:,1);
xfCd6 = xf6(:,3);
xfCl6 = xf6(:,2);
figure
plot(alpha,Cl_6,'g--*',xfalpha6,xfCl6(1:29),'b--+')
title('Cl vs alpha, Re6e5')
xlabel('alpha')
ylabel('Cl')

figure
errorbar(alpha,Cl_6,err6,'r-')
title('Cl vs alpha, Re6e5')
xlabel('alpha')
ylabel('Cl')

figure
count2 = 1;
count3 = -6;
alphaD = [-6:3:24]';
Pdiffsave = zeros(200,11);
while count2<=11
    filechange = '6e5dragdata/alpha_%d.dat';
    pRe6drag{count2}=sortrows(load(sprintf(filechange,count3)));
    count3 = count3+3;
    count = 1;
while count<=200
    p6s = pRe6drag{count2}(count,3:12);
    p6savg = mean(p6s)*47.88;
    p6wake(count,1) = p6savg+patm;
    count = count+1;
end
    qinf6 = (pRe6drag{count2}(1,2).*47.88);
    vinf6 = sqrt(2*qinf6/rho);
    Pfstotal = qinf6+patm;
    Y = pRe6drag{1}(:,1);
    
    Pdiff = (Pfstotal-p6wake)./qinf6;
    Pdiffsave(:,count2) = Pdiff.*qinf6;
    Pdiffsave(Pdiffsave<0)=0;
    
    Cdcorr6(count2) = -trapz(Pdiffsave(:,count2),Y).*1/(c*qinf6);
    errcd6(count2) = std(Cdcorr6)./11;
    plot(Pdiff,Y./c)
        title('Pressure differential across wake, Re6e5')
    xlabel('Pressure differential')
    ylabel('Wake-rake location/c')
    hold on
count2 = count2+1;
end
hold off

figure
errorbar(alphaD,Cdcorr6,errcd6,'r-')
title('Cd vs alpha, Re6e5')
xlabel('alpha')
ylabel('Cl')

figure
plot(alphaD,Cdcorr6,'b--x',xfalpha6,xfCd6(1:29),'k--o')
title('Cd vs alpha, Re6e5')
xlabel('alpha')
ylabel('Cd')

figure
plot(Cdcorr6,Cl_6(1:3:31),'b--x',xfCd6,xfCl6,'k--o')
title('Cd vs Cl, Re6e5')
xlabel('Cd')
ylabel('Cl')
%%
figure
count2 = 1;
count3 = -6;
p7posavg = zeros(44,1);
alpha_7 = [-6:23];
while count2<=29
    filechange = '7e5pressures/alpha_%d.dat';
    pRe7{count2}=load(sprintf(filechange,count3));
    count3 = count3+1;
    count = 1;
while count<=44
    pinf7pos = pRe7{count2}(:,4).*47.88;
    vinf7 = sqrt(2*pinf7pos/rho);
    p7pos = pRe7{count2}(count,5:14);
    p7pos = mean(p7pos)*47.88;
    p7posavg(count) = p7pos;
    count = count+1;
end
Cp7pos = (p7posavg)./(.5*rho*vinf7.^2);
X_7 = pRe7{count2}(:,2);
Y_7 = smooth(Cp7pos);
plot(X_7./c,Y_7,'c')
title('Cp vs x/c, Re7e5')
xlabel('x/c')
ylabel('Cp')
set(gca,'Ydir','reverse');
grid on
hold on
count2 = count2+1;
intCpl7 = Y_7(1:23,1);
intCpu7 = Y_7(24:44,1);
Cl_7(count2) = trapz(X_7,Cp7pos)./c;
err7(count2) = std(Cl_7)./29;
end 

xf7 = importdata('red-airfoil_7e5.pol');
xfalpha7 = xf7(:,1);
xfCd7 = xf7(:,3);
xfCl7 = xf7(:,2);

figure
plot(alpha_7,Cl_7,'g--*',xfalpha7,xfCl7,'b--+')
title('Cl vs alpha, Re7e5')
xlabel('alpha')
ylabel('Cl')

figure
errorbar(alpha_7,Cl_7,err7,'r-')
title('Cl vs alpha, Re7e5')
xlabel('alpha')
ylabel('Cl')

figure
count2 = 1;
count3 = -6;
alpha_7 = [-6:3:21]';
Pdiffsave = zeros(200,10);

while count2<=10
    filechange = '7e5dragdata/alpha_%d.dat';
    pRe7drag{count2}=sortrows(load(sprintf(filechange,count3)));
    count3 = count3+3;
    count = 1;
while count<=200
    p7s = pRe7drag{count2}(count,3:12);
    p7savg = mean(p7s)*47.88;
    p7wake(count,1) = p7savg+patm;
    count = count+1;
end
    qinf7 = (pRe7drag{count2}(1,2).*47.88);
    vinf7 = sqrt(2*qinf7/rho);
    Pfstotal = qinf7+patm;
    Y = pRe7drag{1}(:,1);
    
    Pdiff = (Pfstotal-p7wake)./qinf7;
    Pdiffsave(:,count2) = Pdiff.*qinf7;
    Pdiffsave(Pdiffsave<0)=0;
    
    Cdcorr7(count2) = -trapz(Pdiffsave(:,count2),Y).*1/(c*qinf7);
    errcd7(count2) = std(Cdcorr7)./10;
    
    plot(Pdiff,Y./c)
    title('Pressure differential across wake, Re7e5')
    xlabel('Pressure differential')
    ylabel('Wake-rake location/c')
    hold on
count2 = count2+1;
end
hold off

figure
errorbar(alpha_7,Cdcorr7,errcd7,'r-')
title('Cd vs alpha, Re7e5')
xlabel('alpha')
ylabel('Cl')


figure
plot(alpha_7,Cdcorr7,'b--x',xfalpha7,xfCd7,'k--o')
title('Cd vs alpha, Re7e5')
xlabel('alpha')
ylabel('Cd')

figure
plot(Cdcorr7,Cl_7(1:10),'b--x',xfCd7,xfCl7,'k--o')
title('Cd vs Cl, Re7e5')
xlabel('Cd')
ylabel('Cl')