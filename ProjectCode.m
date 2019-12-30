%%
clc; clear;
gamma = 1.4;
P_atm = 14.696;%psi
%% get data from files
allData = cell(1, 3);
image_details = cell(1, 3);


for b = 1:3
	if b == 1
		%matrix index 1 is calibration
		allData{1, b} = [-10.4, -9.3, -8.2, -6.9, -5.8, -4.8, -4, -3.2, -2.5, -1.9, -1.4, -1, -0.8, -0.5, -0.3, 0]; %gathered data - first value -11.6,  omitted
		image_details{b} = dir('calibration/DSC_*');
	else
		%matrix index 2 is block 0818... and index 3 is 2102
		B = {''; '0818'; '2102'};
		
		filename = sprintf('block%s_theta00/data.dat', B{b, 1});
		allData{1, b} = importdata(filename);
		
		%gets info on the images to import
		image_details{b} = dir(sprintf('block%s_theta00/DSC_*', B{b, 1}));
	end
	
	%create matrix for all images
	for i = 1:length(image_details{b})
		images{b}{i} = rgb2gray(imread([image_details{b}(i).folder, '\', image_details{b}(i).name])); %import as grayscale...
	
		%for calibration
        for L = 1:5
		row_begin = randi([1725,2735]); row_end = row_begin+25; columns_begin = randi([1600,2175]); columns_end = columns_begin+25; %CHANGE THIS DEPENDING ON DATASET
		I_mean{b}(L,i) = mean(mean(images{b}{i}(row_begin:row_end, columns_begin:columns_end)));
        end
        I_mean{b}(i) = mean(I_mean{b}(i));
    end
    
end

%% Image processing


%Step 1 - create calibration curve

%find P_0 = atmospheric pressure
P_0 = P_atm * 6894.757;%psi to Pa

I_0 = I_mean{1}(1,14); %mean Intensity value at the image where gauge pressure = 0, meaning pressure = P_0 = P_atm

%convert P values in allData to be absolute and SI
P{1} = (allData{1} + P_atm) * 6894.757;

%find ratio P/P_0 for every value
P_ratio{1} = P{1} / P_0;

%intensity ratio  =  Imean (for each image) / I_0 (I when P = P_0)
I_ratio{1} = I_mean{1} / I_0;
figure
for i = 1:10
for L=1:5
%get straight line of best fit
 c(L,1:2) = polyfit(P_ratio{1}(1:14), I_ratio{1}(L,:), 1);
 l = polyval(c(L,:), P_ratio{1});

%plot calibration curve
hold on
plot(P_ratio{1}(1:14), I_ratio{1}(L,i),'*');
plot(P_ratio{1}, l);
end
end


%plot(P_ratio{1}, l);



hold off;
%NOTE: right now we only are using one value from each image, i think in the example he uses 10















