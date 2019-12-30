%  time (s) | P_01 (psi) | P_1 (psi) | P_02 (psi) | Patm (psi) | T_01 (degree F)
% Column definitions in 'manufacturer_calibration_curve_supersonic_tunnel.txt': Block Setting | Mach Number
% - The file names correspond to the block settings at which the data was collected.
% - All recorded pressures are gauge pressures. Remember to add the atmospheric pressure during your analysis.
% - Only consider data beyond 2.5 seconds of run time. This is due to the fact that it takes the flow that amount of time 
%   to reach the pre-set freestream stagnation pressure (P_01) of 60 psi.
% - Suggestion - In order to get the solution for the Mach number from the Rayleigh Pitot equation, you can use, fsolve(), 
%   the in-built Newton-iterator in MATLAB. More information on the implementation of the same can be found on MATLAB's on-line documentation.
% - The report needs to be in the AIAA prescribed format.
% - All results in your final report need to be converted to SI units.
% - As discussed in the lab, focus on critial analysis. Point out and discuss interesting trends (if any) in the data.
% - In order to get the extra credit, you will have to present detailed analysis and plots. 
fileref = 400;
gam = 1.4;
blockN = zeros(12,1);
for A=1:12  
    blockN(A) = fileref;
    filenameformat = 'calibration_%d.txt';
    filename = sprintf(filenameformat,fileref);
    if A==1
    Data = dlmread(filename,'\t',127,0);
    else
        newfile = dlmread(filename,'\t',127,0);
        Data = [Data; newfile];
    end
    fileref = fileref + 200;
end
Patm = mean(Data(:,5));
AdjPressures = (Data(:,2:4)+Patm)*6894.76; %psi to pa
avPre = zeros(12,3);
refL = 1;
refU = 126;
for b=1:12
avPre(b,1) = mean(AdjPressures(refL:refU,1));
avPre(b,2) = mean(AdjPressures(refL:refU,2));
avPre(b,3) = mean(AdjPressures(refL:refU,3));
refL = refL+126;
refU = refU+126;    
end

M_ise = sqrt(2*((avPre(:,1)./avPre(:,2)).^((gam-1)/gam)-1)/(gam-1));
syms M
Mach = zeros(12,1);

for C=2:12
eqn = (((gam+1)^2*M^2)/(4*gam*M^2-2*(gam-1)))^(gam/(gam-1))*((1-gam+2*gam*M^2)/(gam+1))== avPre(C,3)./avPre(C,2);
Mach(C) = double(vpasolve(eqn,M,2.3));
end

figure
plot(blockN(3:12,1),Mach(3:12,1))
xlabel('Block #') 
ylabel('Mach #')
title('Rayleigh-Pitot relation')

figure
plot(blockN(3:12,1),M_ise(3:12,1))
xlabel('Block #')
ylabel('Mach #')
title('Isentropic relation')

figure
plot(blockN(1:12,1),Mach(1:12,1))
xlabel('Block #') 
ylabel('Mach #')
title('Rayleigh-Pitot relation')

figure
plot(blockN(1:12,1),M_ise(1:12,1))
xlabel('Block #')
ylabel('Mach #')
title('Isentropic relation')
