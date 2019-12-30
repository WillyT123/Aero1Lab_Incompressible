
%|  Time(s)  | P01  |  P1  |  Patm  |  T01  |
tic
clc
close all
imtextref = [6 16];
cd all-data
folders = dir(pwd);
folders(1:2) = [];
gam=1.4;
avgData = zeros(35,5); alpha = zeros(1,35); angleBeta = zeros(35,5);
bn = avgData; X=zeros(537,1); Y=X; convert=zeros(251,5);

for n=1:length(folders)
    bn(n)=str2double(folders(n).name(6:9));
    alpha(n)=str2double(folders(n).name(16:17));
    cd(folders(n).name);
    files = dir(pwd);
    files(1:2)=[];
    fileopen = fopen('data.dat','r');
    convert = cell2mat(textscan(fileopen,repmat('%f',1,5),'headerlines',2));
    datatrim = find(2.5<convert(:,1),1);
    avgData(n,:)=mean(convert(datatrim:end,:),1);
    fclose('all');   
    for x=1:5  
        pic=imread(files(x+8).name);
        variation1 = diff(pic);
        variation2 = diff(pic, 2);
        edgeTol1 = max(mean(variation1));
        edgeTol2 = max(mean(variation2));
        threshold = mean([edgeTol1, edgeTol2])*0.175;
        edges = edge(pic, 'canny', threshold);
        logical = edges & cumsum(edges, 1) == 1;
        [Y, X] = find(logical);
        poly = polyfit(X, Y, 1); % X and Y are flipped
        angleBeta(n,x) = atand(poly(1));
  
    end
    cd ../
end
cd ../

avAngles=zeros(35,1);
for J=1:35
avAngles(J)=mean(angleBeta(J,:));     %Betas
end

%%
thetas = linspace(0,20,5);
bta = zeros(7,5);
t = 1;
r = 5;
for I=1:7
    bta(I,:) = avAngles(t:r);
    t=t+5;
    r=r+5;
 end
    

Patm = mean(avgData(1,4));
AdjP = (avgData(:,2:3)+Patm)*6894.76;
M_ise = sqrt((nthroot((AdjP(:,1)./AdjP(:,2)),(gam/(gam-1)))-1)/((gam-1)/2));

Pmach = zeros(7,1);
t = 1;
r = 5;
for I=1:7
    Pmach(I) = mean(M_ise(t:r,1));
    t=t+5;
    r=r+5;
end

bn = [818, 895, 1012, 1286, 1667, 2102, 2465];
M_i = 1.82e-7*bn.^2-1.3e-3*bn+3.90;

M = zeros(7,5);
for i = 1:7
    for j=1:5
        M(i,j) = (sind((bta(i,j))).^2-(1.2./(1+cotd(thetas(j)).*cotd(bta(i,j))))).^-0.5;
    end
end

m = zeros(1,7);
for k=1:7
    if k==6
    m(k) = mean(M(k,2:5));
    else
    m(k) = mean(M(k,:));
    end
    
end
figure
plot(bn,m,'c--')
hold on
plot(bn,M_i,'b--')
hold on
plot(bn,Pmach,'g--')
hold on
title('Block # vs Mach #')
xlabel('Block #')
ylabel('Mach #')
grid on
%%

figure
for I=1:7
    if I==6
        plot(thetas(2:5),bta(I,2:5),'g')
        hold on
    else
    plot(thetas,bta(I,:),'g') 
    hold on
    end
end


syms B
Beta = zeros(7,5);
thetas(1) = 1e-20;
thetas(5) = 19.99999999999999;
bn = [818, 895, 1012, 1286, 1667, 2102, 2465];
M_i = 1.82e-7*bn.^2-1.3e-3*bn+3.90;
for c = 1:7
for C=1:5
eqn = tand(B-thetas(C))/tand(B) == (2+(gam-1)*M_i(c)^2*sind(B)^2)/((gam+1)*M_i(c)^2*sind(B)^2);
Beta(c,C) = double(vpasolve(eqn,B,35));
end
if c==7
    plot(thetas(1:4),Beta(c,1:4),'r') 
else
plot(thetas,Beta(c,:),'r')
title('Beta vs Theta')
xlabel('Theta')
ylabel('Beta')

grid on
hold on
end
end
toc