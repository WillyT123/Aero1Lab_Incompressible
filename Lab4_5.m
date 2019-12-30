gam  = 1.4;
Ai = 3.17e-5; %m^2
To = 22+273; %deg K


%Lab 4
% Tap Number | Tap Axial Positon (inches) | Nozzle Area Ratio (A/Ai) | P_static (psi) | P_O (psi) | Mass Flow Rate (slugs/second) | P_atm (psi) 

%Lab 5

for A=1:24  

    filenameformat = 'run_%d.txt';
    filename = sprintf(filenameformat,A);
    if A==1
    Data4 = dlmread(filename,'\t');
    else
        newfile = dlmread(filename,'\t');
        Data4 = [Data4; newfile];
    end
 
        filenameformat = 'test_%d.txt';
    filename = sprintf(filenameformat,A);
    if A==1
    Data5 = dlmread(filename,'\t');
        
    else
        newfile = dlmread(filename,'\t');

        Data5 = [Data5; newfile];
    end
    
end
    
    
   Data4(:,1) = [];
   Data5(:,1) = [];
   
   
for L=1:240
    
   if Data5(L,1)==0 
    Data5(L,:)=[];
   end
    
end 


   Data4(:,2) = Data4(:,2)*.0254;
   Data5(:,2) = Data5(:,2)*.0254;
   
   Data4(:,4:5) = Data4(:,4:5)*6894.76+101325;
   Data5(:,4:5) = Data5(:,4:5)*6894.76+101325;            %Converting from imperial to metric units and accounting for gauge pressure measurements
   Data4(:,7) = Data4(:,7)*6894.76+101325;
   Data5(:,7) = Data5(:,7)*6894.76+101325;
   
   Data4(:,6) = Data4(:,6)*14.5939;
   Data5(:,6) = Data5(:,6)*14.5939;
   
   %%
   
   
   Prat4 = zeros(240,1);
   Prat5 = zeros(240,1);
   Pb_rat4 = zeros(24,1);
   Pb_rat5 = zeros(24,1);
   Pt_rat5 = zeros(24,1);
   Pe_rat5 = zeros(24,1);
   MFP_4 = zeros(24,1);
   MFP_5 = zeros(24,1);
  
   
   for B=1:240 
    Prat4(B,1) = Data4(B,4)/Data4(B,5);
    Prat5(B,1) = Data5(B,4)/Data5(B,5);
   end
    
 for C=1:24
     Pb_rat4(C,1) = Prat4(C*10);
     Pb_rat5(C,1) = Prat5(C*10);
     Pt_rat5(C,1) = Prat5(C*10-7);
     Pe_rat5(C,1) = Prat5(C*10-1);
     MFP_4(C,1) = Data4(C*10,6)*sqrt(To)/(Data4(C*10-1,3)*Ai*Data4(C*10,5));
     MFP_5(C,1) = Data5(C*10,6)*sqrt(To)/(Data5(C*10-7,3)*Ai*Data5(C*10,5));
 end

Mach4 = zeros(240,1);
Mach5 = zeros(240,1);

for C=1:240

Mach4(C,1) = abs(sqrt(2/(gam-1)*(Prat4(C)^((gam-1)/-gam)-1)));
Mach5(C,1) = abs(sqrt(2/(gam-1)*(Prat5(C)^((gam-1)/-gam)-1)));

end

figure
for D=1:24
    hold on
    Lref = D*10-9;
    Uref = D*10;
    
plot(Data4(Lref:Uref,2),Mach4(Lref:Uref),'g')
xlabel('Tap dist. (m)')
ylabel('Mach #')
title('Converging Nozzle, Mach Numbers v Distance')
end
plot(linspace(0,.18,10),ones(10),'k--')
    hold off
%%
line = ones(10,1)*0.5283;
figure
for D=1:24
    hold on
    Lref = D*10-9;
    Uref = D*10;
    
plot(Data5(Lref:Uref,2),Mach5(Lref:Uref),'m')
ylabel('Mach #')
title('Converging-Diverging Nozzle, Mach Numbers v Distance')
end
plot(linspace(0,0.25,10),ones(10),'k--')
    hold off

figure
for E=1:24
    hold on
    Lref = E*10-9;
    Uref = E*10;

plot(Data5(Lref:Uref,2),Prat5(Lref:Uref),'b')
xlabel('Tap dist. (m)')
ylabel('P/Po')
title('Converging-Diverging Nozzle, P/Po v Distance')
end
plot(linspace(0,.25,10),line,'k--')
    hold off
%%
figure
for F=1:24
    hold on
    Lref = F*10-9;
    Uref = F*10;
    
plot(Data4(Lref:Uref,2),Prat4(Lref:Uref),'r')
xlabel('Tap dist. (m)')
ylabel('P/Po')
title('Converging Nozzle, P/Po v Distance')
end
plot(linspace(0,.18,10),line,'k--')
hold off
%%
line2 = ones(10,1)*0.5283;
figure
hold on
plot(Pb_rat4,MFP_4,'c*')
plot(line2,linspace(0,0.25,10),'k--')
xlabel('Pb/Po')
ylabel('Mass Flow Parameter')
title('Converging Nozzle, MFP v Back Pressure Ratio')
hold off

figure
hold on
plot(Pb_rat5,MFP_5,'bo')
plot(line2,linspace(0,0.25,10),'k--')
xlabel('Pb/Po')
ylabel('Mass Flow Parameter')
title('Converging-Diverging Nozzle, MFP v Back Pressure Ratio')
hold off
%%
figure
hold on
plot(Pb_rat5,Pe_rat5,'g+')
plot(line2,linspace(0,1,10),'k--')
plot(linspace(0.2,1,10),line2,'k--')
xlabel('Pb/Po')
ylabel('Pe/Po')
title('Converging-Diverging Nozzle, Pe/Po v Pb/Po')
hold off

figure
hold on
plot(Pb_rat5,Pt_rat5,'md')
plot(line2,linspace(0.5,1,10),'k--')
plot(linspace(0.2,1,10),line2,'k--')
xlabel('Pb/Po')
ylabel('Pt/Po')
title('Converging-Diverging Nozzle, Pt/Po v Pb/Po')