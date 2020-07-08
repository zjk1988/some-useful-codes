clear all;
clc;
addpath('E:\experiments\CT_reconstruction\bin');
imgN = 256;
dImg = 0.6641 / 100;
detN = 512;
detL = 0.7200*512 / 100;
angN = 1024;
scanR = 500 / 100 / 2;
detR = 500 / 100 / 2;
% auto calculation
dAng = 0.006134*180/pi;
lowdose = 2.5e4;
lowdose_file = '0.025';
dDet = detL / detN;
dedge = detL / 2 -dDet / 2;
dd = dedge * scanR / sqrt(dedge^2 + (scanR+detR)^2);
mask1 = ((-imgN/2 + 1/2):1:(imgN/2 - 1/2)).^2;
mask = zeros(-imgN,-imgN);
for iii=1:imgN
    mask(:,iii) = sqrt(mask1 + mask1(iii)) * dImg;
end
mask = (mask <= dd);

vrange = (1000 + 1000) / 1000 * 1.92;
show_win = [(-160 + 1000) / 1000 * 1.92 /vrange, (240 + 1000) / 1000 * 1.92 /vrange];
% 
detR1 = 0;

dDet1 = dDet * scanR/(scanR+detR);
detL1 = detL * scanR/(scanR+detR);

% conv kernel R-L
i=1;
for n=-detN+1:detN-1
    if(mod(n,2)==1)
        h(i)=-1/(pi*pi.*n.*n*dDet1^2);
    else
        if(n==0)
            h(i)=1/(4*dDet1^2);
        else
            h(i)=0;
        end
    end
    i=i+1;
end
% 
% % pre-weighting
w = ((-detN/2+1/2)*dDet1:dDet1:(detN/2-1/2)*dDet1)';
W = scanR./sqrt(scanR^2+w.^2);

label_dir = strcat('train\label_single');
noisy_dir = strcat('train\projection_',lowdose_file);
input_dir = strcat('train\input_',lowdose_file);
img_dir = strcat('train\img_',lowdose_file);
if ~exist(label_dir)
    mkdir(label_dir);
end
if ~exist(noisy_dir)
    mkdir(noisy_dir);
end
if ~exist(input_dir)
    mkdir(input_dir);
end
if ~exist(img_dir)
    mkdir(img_dir);
end

for ii = 1
    path = fullfile('train','noisy_free');
    img_list = dir(path);
    for tt = 1:length(img_list)
        if strcmp(img_list(tt).name, '.') || strcmp(img_list(tt).name, '..')
            continue;
        end
        file_name = fullfile(path,img_list(tt).name);
        orig_name = strrep(file_name,'noisy_free','label');
        load(orig_name, 'data');
        orig = data;
        orig = single(orig / vrange);
        load(file_name, 'data');
        P_FBP = data';
        %convolution
        P_FBP = addPoissNoise(P_FBP, lowdose, 10);
        P_FBP = P_FBP / vrange;
        for i=1:size(P_FBP,2)
            G_FBP(:,i) = conv(dDet1*W.*P_FBP(:,i),h,'same');
        end
        I_FBP = xrec(G_FBP,'E:\experiments\CT_reconstruction\demos\default.cfg',...
            'Scanning_Type','FAN_ED',...
            'XRay_To_Rotation_Center',scanR,'Detector_To_Rotation_Center',detR1,...
            'Total_Num_Of_Angles',angN,'First_Angle_Degree',0,'Angle_Increment',dAng,...
            'Length_Of_Detector',detL1,'Num_Of_Detector_Units',detN,...
            'Rows_Of_Image',imgN,'Columns_Of_Image',imgN,...
            'Pixel_Width',dImg,'Pixel_Height',dImg,...
            'Reconstruction_Algorithm','BACKPRJ',...
            'Regularization',0);
        I_FBP = I_FBP / 2 .* mask;
        figure(2);imshow(I_FBP, show_win);title('FBP Reconstuction');
        label_name = fullfile(label_dir, img_list(tt).name);
        input_name = fullfile(input_dir, img_list(tt).name);
        prj_name = fullfile(noisy_dir, img_list(tt).name);
        img_name = fullfile(img_dir, strcat(img_list(tt).name(1:end-4),'.jpg'));
        data = orig;
        save(label_name,'data');
        data = single(I_FBP);
        save(input_name,'data');
        data = single(P_FBP');
        save(prj_name,'data');
        I_FBP = uint8((I_FBP-0.42)/(0.62-0.42) * 255);
        imwrite(I_FBP, img_name);
    end
end

label_dir = strcat('test\label_single');
noisy_dir = strcat('test\projection_',lowdose_file);
input_dir = strcat('test\input_',lowdose_file);
img_dir = strcat('test\img_',lowdose_file);
if ~exist(label_dir)
    mkdir(label_dir);
end
if ~exist(noisy_dir)
    mkdir(noisy_dir);
end
if ~exist(input_dir)
    mkdir(input_dir);
end
if ~exist(img_dir)
    mkdir(img_dir);
end

for ii = 1
    path = fullfile('test','noisy_free');
    img_list = dir(path);
    for tt = 1:length(img_list)
        if strcmp(img_list(tt).name, '.') || strcmp(img_list(tt).name, '..')
            continue;
        end
        file_name = fullfile(path,img_list(tt).name);
        orig_name = strrep(file_name,'noisy_free','label');
        load(orig_name, 'data');
        orig = data;
        orig = single(orig / vrange);
        load(file_name, 'data');
        P_FBP = data';
        %convolution
        P_FBP = addPoissNoise(P_FBP, lowdose, 10);
        P_FBP = P_FBP / vrange;
        for i=1:size(P_FBP,2)
            G_FBP(:,i) = conv(dDet1*W.*P_FBP(:,i),h,'same');
        end
        I_FBP = xrec(G_FBP,'E:\experiments\CT_reconstruction\demos\default.cfg',...
            'Scanning_Type','FAN_ED',...
            'XRay_To_Rotation_Center',scanR,'Detector_To_Rotation_Center',detR1,...
            'Total_Num_Of_Angles',angN,'First_Angle_Degree',0,'Angle_Increment',dAng,...
            'Length_Of_Detector',detL1,'Num_Of_Detector_Units',detN,...
            'Rows_Of_Image',imgN,'Columns_Of_Image',imgN,...
            'Pixel_Width',dImg,'Pixel_Height',dImg,...
            'Reconstruction_Algorithm','BACKPRJ',...
            'Regularization',0);
        I_FBP = I_FBP / 2 .* mask;
        figure(2);imshow(I_FBP, show_win);title('FBP Reconstuction');
        label_name = fullfile(label_dir, img_list(tt).name);
        input_name = fullfile(input_dir, img_list(tt).name);
        prj_name = fullfile(noisy_dir, img_list(tt).name);
        img_name = fullfile(img_dir, strcat(img_list(tt).name(1:end-4),'.jpg'));
        data = orig;
        save(label_name,'data');
        data = single(I_FBP);
        save(input_name,'data');
        data = single(P_FBP');
        save(prj_name,'data');
        I_FBP = uint8((I_FBP-0.42)/(0.62-0.42) * 255);
        imwrite(I_FBP, img_name);
    end
end