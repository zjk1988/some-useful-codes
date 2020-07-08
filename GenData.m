% clc; clear; close all;
% function generate_simulated_suitcases_Amat_gt(start)
format long;
rng(512,'twister'); % fixing seed
start = 1
% if ~exist('A','var')
%     load('Amatrix.mat');
% end
%% loading LACs and spectral response
[num, str, raw] = xlsread('LAC_final_scase.xlsx');
energy = num(:,1)*1e3;
% converting LACs from /cm to /mm
num(:,2:end) = num(:,2:end)/10*475/512;
nE = 131;%size(energy,1);
data1 = 'concat_suitcases_gt_Amat_512/';
mkdir(data1);

% Spectral response
W = zeros(nE,2);
M = dlmread('imatron_spectrum/imatron_spectrum_95kV.txt');
W(1:86,1) = M(:,2);
M = dlmread('imatron_spectrum/imatron_spectrum_130kV.txt');
W(1:121,2) = M(:,2);
%% simulation parameters
numTrain = 10000;
numTest = 100;
numMetal = 83;
numNonMetal = 10;
maxNonMetal = 50;

num_photons_130=170000;
num_photons_95=180000;
angles = 0:0.25:179.9;
imSize = 512;

% astra
vol_geom = astra_create_vol_geom(imSize, imSize);
proj_geom = astra_create_proj_geom('parallel', 0.5, 1024, linspace2(0,pi,720));

% computing LACs for air
load('elements_info.mat');
load('compounds_xls_table.mat');
% air dry, near sea level
air = 1/10*1.205e-03 *( elements_info(6).mass_attenuation*0.000124 + elements_info(7).mass_attenuation*0.755268 + elements_info(8).mass_attenuation*0.231781  + elements_info(18).mass_attenuation*0.012827);
% air equivalent plastic
% air = 1/10* 1.760 *( elements_info(1).mass_attenuation*0.024681 + elements_info(6).mass_attenuation*0.501610 + elements_info(8).mass_attenuation*0.004527  + elements_info(9).mass_attenuation*0.465209 + elements_info(14).mass_attenuation*0.003973);
% neoprene to simulate clothes
neoprene = num(:,11);
%% generating simulated data

% randomly select number of non-metal items to be placed in the phantom in the range [0,3]
num_items = randi([0 maxNonMetal],numTrain+numTest,1);
load('scase_bg1.mat');
temp = zeros(size(scase_bg));
for i = 1:nE
    scase = scase_bg(:,:,i);
    scase(136:378,14:498) = neoprene(i)/10;
    idcs = find(scase==0);
    scase(idcs) = air(i);
    temp(:,:,i) = scase;
end
scase_bg = temp;

for n = 1:numTrain
    n
    tic
    maxMetal = 5;
    % randomly select non-metal items
    nonMetal = randi([1 numNonMetal],num_items(n),1);
    % randomly select nmetal items
    metal = randi([1 numMetal],maxMetal,1);
    % randomly select shape of jth item: 1 rectangle, 2 ellipse
    shape_type = randi([1 2],num_items(n)+maxMetal,1);
    % randomy select min and max axis length of rectangle/ellipse
    axis_len = randi([5 150],num_items(n),2);
    axis_len_m = randi([5 100],maxMetal,2);
    axis_len = [axis_len;axis_len_m];
    % randomly select center points of ellipses/rectangles
    a1 = rand_minsep(num_items(n)+maxMetal, 40, 512-40, 10*ones(num_items(n)+maxMetal,1)); % <- here is your random array of n points for cols
    %update number of objects that were actually placed in the suitcase
    num_items(n) = length(a1) - maxMetal;
    
    a2 = rand_minsep(num_items(n)+maxMetal, 40, 512-40, 10*ones(num_items(n)+maxMetal,1)); % <- here is your random array of n points for rows
    swap = max(0,min(length(a1),length(a2)) - maxMetal);
    maxMetal = min(maxMetal, min(length(a1),length(a2)) );
    num_items(n) = swap;
    centers = zeros(num_items(n)+maxMetal,2);
    for i=1:num_items(n)+maxMetal
        centers(i,:) = [a1(i), a2(i)];
        if (ceil(axis_len(i,1)/2) + centers(i,1)) > imSize
            axis_len(i,1) = imSize - centers(i,1) - 1;
        elseif (centers(i,1) - ceil(axis_len(i,1)/2)) < 1
            axis_len(i,1) = centers(i,1) - 1;
        end
        
        if (ceil(axis_len(i,2)/2) + centers(i,2)) > imSize
            axis_len(i,2) = imSize - centers(i,2) - 1;
        elseif (centers(i,2) - ceil(axis_len(i,2)/2)) < 1
            axis_len(i,2) = centers(i,2) - 1;
        end
        
    end
    n_name = sprintf('%s/sino95_n_%d_m5_m0_num_%d.png',data1,n,num_items(n));
%     if exist(n_name,'file') == 2
%         continue;
%     end
    
    % to enable multiple parallel runs using qsub
    if n < 1
        continue;
    elseif n > 1000
        break;
    end
    
    % genereating phantom
    
    phantom = scase_bg;%zeros(imSize,imSize,nE);
    for j = 1:num_items(n)
        if shape_type(j) == 2
            r = ceil(min(axis_len(j,:)/2));
            for i = 1:nE
                phantom(:, :,i) = midpointCircle(phantom(:, :,i), r , centers(j,1), centers(j,2), num(i, nonMetal(j)+1));
            end
            %             rec_gt = midpointCircle(rec_gt, r , centers(j,1), centers(j,2), num(51, nonMetal(j)+1));
        elseif shape_type(j) == 1
            
            for i = 1:nE
                phantom(max(1,centers(j,2) - ceil(axis_len(j,2)/2)): min(imSize,centers(j,2) + ceil(axis_len(j,2)/2)), max(1,centers(j,1) - ceil(axis_len(j,1)/2)): min(imSize,centers(j,1) + ceil(axis_len(j,1)/2)),i) = num(i, nonMetal(j)+1);
            end
            %             rec_gt(centers(j,2) - ceil(axis_len(j,2)/2): centers(j,2) + ceil(axis_len(j,2)/2), centers(j,1) - ceil(axis_len(j,1)/2): centers(j,1) + ceil(axis_len(j,1)/2)) = num(51, nonMetal(j)+1);
        end
    end
    
    for k = 1:maxMetal
        if shape_type(num_items(n)+k) == 2
            r = ceil(min(axis_len(size(axis_len,1)-maxMetal+k,:)/2));
            for i = 1:nE
                phantom(:, :,i) = midpointCircle(phantom(:, :,i), r , centers(num_items(n)+k,1), centers(num_items(n)+k,2), air(i));
            end
        elseif shape_type(num_items(n)+k) == 1
            for i = 1:nE
                phantom(max(1,centers(num_items(n)+k,2) - ceil(axis_len(size(axis_len,1)-maxMetal+k,2)/2)): min(imSize,centers(num_items(n)+k,2) + ceil(axis_len(size(axis_len,1)-maxMetal+k,2)/2)), max(1,centers(num_items(n)+k,1) - ceil(axis_len(size(axis_len,1)-maxMetal+k,1)/2)): min(imSize,centers(num_items(n)+k,1) + ceil(axis_len(size(axis_len,1)-maxMetal+k,1)/2)),i) = air(i);
            end
        end
    end
    
    rec_gt = phantom(:,:,51);
    
    proj = zeros(1024,720, 121);
    for i = 1:121
        %             tic;
        slice = phantom(:,:,i);
        %             temp = A*slice(:);
        [sinogram_id, sinogram] = astra_create_sino_gpu(slice*2, proj_geom, vol_geom);
        proj(:,:,i) = sinogram';
        astra_mex_data2d('delete', sinogram_id);
        %             toc
    end
    
    % create polyenergetic sinograms
    weights = energy*0.0025954;
    [sg_130, sg_noise_130, sg_ln_130, sg_ln_noise_130]=combinePolyEDataENoise(proj(:,:,1:121), energy(1:121), W(1:121,2), weights(1:121) , num_photons_130);
    %         [sg_95, sg_noise_95, sg_ln_95, sg_ln_noise_95]=combinePolyEDataENoisePoisson(proj(:,:,1:86), energy(1:86), W(1:86,2), weights(1:86), num_photons_95,u);
    m0 = sg_ln_noise_130*4.2296;
    %         rec_gt = gather(rec_gt);
    
    %         phName = sprintf('%s/ph_n_%d_m_%d_num_%d.tiff',data1,n,m,num_items(n));
    %         imwrite(rec_gt,phName);
    %         phName = sprintf('%s/ph_n_%d_m_%d_num_%d.mat',data1,n,m,num_items(n));
    %         save(phName,'rec_gt');
    %% metal related part of the code
    for m = 1:maxMetal
        s95Name = sprintf('%s/sino95_n_%d_m%d_m0_num_%d.png',data1,n,m,num_items(n));
        %         if exist(s95Name,'file') == 2
        %         	continue;
        %         end
        metal = zeros(imSize);
        for k = 1:m
            if shape_type(num_items(n)+k) == 2
                r = ceil(min(axis_len(size(axis_len,1)-maxMetal+k,:)/2));
                metal = midpointCircle(metal, r , centers(num_items(n)+k,1), centers(num_items(n)+k,2), 1);
            elseif shape_type(num_items(n)+k) == 1
                metal(max(1,centers(num_items(n)+k,2) - ceil(axis_len(size(axis_len,1)-maxMetal+k,2)/2)): min(imSize,centers(num_items(n)+k,2) + ceil(axis_len(size(axis_len,1)-maxMetal+k,2)/2)), max(1,centers(num_items(n)+k,1) - ceil(axis_len(size(axis_len,1)-maxMetal+k,1)/2)): min(imSize,centers(num_items(n)+k,1) + ceil(axis_len(size(axis_len,1)-maxMetal+k,1)/2))) = 1;
            end
        end
        %         proj = A*metal(:);
        %         proj = reshape(proj, [1024 720]);
        [sinogram_id, sinogram] = astra_create_sino_gpu(metal, proj_geom, vol_geom);
        proj = sinogram';
        astra_mex_data2d('delete', sinogram_id);
        % getting rid of the metal related parts in the projection data
        m1 = m0;
        idcs = find(proj);
        m1(idcs) = 0;
        
        mask = (proj>0);
        m1_m0 = [m0,mask,m1];
        
        %         figure;imagesc(m1_m0);
        m1_m0 = gather(m1_m0);
        imwrite(m1_m0/50,s95Name);
        %         save( s95Name, 'm1_m0', '-v6');
    end
    toc
end