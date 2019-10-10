% Creates .mat file for personal tests on hadamard code, using generated
% Hadamard optical sectioning microscopy patterns in complement interleaved
% fashion.

% GENERATE HADAMARD CODES, COMPLEMENT INTERLEAVED

% Define number of orthogonal locations.
% Usually as n = m - 1 with m a multiple of 4,
% in any case, m>n patterns are generated, with m as low as possible.
% The offset defines how orthogonal locations distribute in space; they
% should be equidistant from each other to minimize scattering crosstalk
nlocations_and_offset = [19 5]; % [n offset]

% optional parameter binningto group projection elements
binning = [1 1];

% Generate a matrix of size [1024 768]. randomization is always the same.
patterns_logical = alp_btd_to_logical(hadamard_patterns_scramble_nopermutation(nlocations_and_offset,binning));

% Interleave with complement in 3rd dimension, store as vectorized movie
toint = {... % to interleave matrices in 3d, place them in a cell to then use the permute-cell2mat-permute-reshape method
     patterns_logical ,...
    ~patterns_logical};
hadamard_patterns = vm(reshape(permute(cell2mat(permute(toint,[1 4 3 2])),[1 2 4 3]),size(toint{1},1),size(toint{1},2),[]));
without_complement = vm(patterns_logical);

% figure(5)
% moviesc(hadamard_patterns)

% Crop the matrix around ROI for easier visualization
hadamard_patterns = hadamard_patterns(1:44,1:80,:);
without_complement = without_complement(1:44,1:80,:);

% Display patterns as imagej stack style. note there are 2*m frames
figure(1)
moviesc(hadamard_patterns)
title 'illumination patterns'

figure(2)
moviesc(without_complement)
title 'illumination patterns without interleaved complements'

% GENERATE BASE IMAGE
% [dim1, dim2, dim3] = size(hadamard_patterns);
% 
% [columnsInImage, rowsInImage] = meshgrid(1:dim2, 1:dim1);
% centerX = 40;
% centerY = 22;
% radius = 10;
% circlePixels = (rowsInImage - centerY).^2 + (columnsInImage - centerX).^2 <= radius.^2;
% % Euclidean Distance Transform
% distanceImage = bwdist(~circlePixels);
% % Normalization
% distanceImage = distanceImage / max(distanceImage(:));
% distanceImage = double(distanceImage);
% 
% figure(3);
% d = vm(distanceImage);
% moviesc(d);

% GENERATE "MOVING" BASE IMAGES
[dim1, dim2, dim3] = size(hadamard_patterns);

timestep1 = generateBaseImageGradientCircle(20,22,10,dim1,dim2);
timestep2 = generateBaseImageGradientCircle(40,22,10,dim1,dim2);
timestep3 = generateBaseImageGradientCircle(60,22,10,dim1,dim2);

distanceImage = cat(3, timestep1, timestep2, timestep3);

figure(3);
d = vm(distanceImage);
moviesc(d);

% ADD NOISE TO EACH FRAME, APPLY HADAMARD PATTERN
reps = 40;
[imgdim1, imgdim2, imgdim3] = size(distanceImage);
image_patterns = zeros(dim1, dim2, dim3*reps*imgdim3);
for k = 1:imgdim3
    for i = (dim3*reps*(k-1))+1:(dim3*reps*k)
        j = mod(i, dim3);
        if j == 0
            j = dim3;
        end
        pattern = hadamard_patterns.data(:,:,j);
        new_img = distanceImage(:, :, k) .* pattern;
        new_img_with_gauss = imnoise(new_img,'gaussian');
        new_img_with_poiss = imnoise(new_img, 'poisson');
        image_patterns(:, :, i) = new_img_with_gauss + new_img_with_poiss - new_img;
    %     disp(mean(image_patterns(12:32, 30:50, i))/std(image_patterns(12:32, 30:50, i)));
    end
end

figure(4);
image_patterns = vm(image_patterns);
moviesc(image_patterns);

hadamard_patterns = double(hadamard_patterns);



% RECONSTRUCT IMAGE
rootdir = fullfile('Examples/Raw data');
%% Code parameters
il = 20*2; % interleaved length containing 2 series of length 20
hl = il/2; % hadamard length, equal to half the interleaved length
hadtraces = hadamard_bincode_nopermutation(hl-1)'*2-1;

%% Format calibration data
horiz_dim = 80;
vert_dim = 44;
%%
    res_fpath = fullfile(rootdir,'results.mat');

    %%
    ncomps = 5;
    rng(0,'twister') % ensures multiple runs of eigs() will return repeatable results
%     resid_mov = dyn_had_mod(image_patterns(:,:), without_complement(:,:), ncomps);
%     resid_mov = vm(reshape(resid_mov.data, vert_dim, horiz_dim, dim3*reps));
%     figure(50);
%     moviesc(resid_mov);
%     [sv, uhad, uref, uu] = dyn_had(resid_mov(:,:), without_complement(:,:), ncomps);
    [sv, uhad, uref, uu] = dyn_had(image_patterns(:,:), without_complement(:,:), ncomps);
    %%
    save(res_fpath,...
        'sv', ...
        'uhad', ...
        'uref', ...
        'uu', ...
        'ncomps')

%% Analysis
% manually defined ROIs to integrate cell excluding overlapping area
rois = {[
   11.5251   27.6025
   18.6321   32.8872
   24.2813   26.6913
   21.5478   18.4909
   25.3747    8.6503
   33.9396    6.8280
   39.4066    1.3610
    4.2358    1.1788
    1.8667    1.3610
    2.2312   23.4112
   11.5251   27.6025
    ],[
   20.0900   43.4567
   20.6367   35.2563
   25.1925   26.8736
   32.6640   24.8690
   43.9624   21.4066
   47.9715   12.6595
   48.1538    5.9169
   53.2563    7.1925
   60.7278   20.3132
   62.5501   32.1583
   61.6390   42.7278
   20.0900   43.4567
    ]};
%% 

% reload reconstructed data
res_fpath = fullfile(rootdir,'results.mat');
load(res_fpath,...
    'sv', ...
    'uhad', ...
    'uref', ...
    'uu', ...
    'ncomps')
% extract time-integrated images 
w50img = mean(vm(uref*sv',image_patterns.imsz)); % widefield image
h50img = mean(vm(uhad*sv',image_patterns.imsz)); % hadamard image
% h50img = 1-mean(vm(uhad*sv',image_patterns.imsz)); % hadamard image

% reload reconstructed data from depth 60 um
res_fpath = fullfile(rootdir,'results.mat');
load(res_fpath,...
    'sv', ...
    'uhad', ...
    'uref', ...
    'uu', ...
    'ncomps')
% extract time-integrated images 
w60img = mean(vm(uref*sv',image_patterns.imsz));
h60img = mean(vm(uhad*sv',image_patterns.imsz));

figure(5);
moviesc(vm(h50img));

% mse_19_5_400_svdtwice = immse(distanceImage, h50img);
% name = 'mse_19_5_400_svdtwice';
% save('Simulations-lower_intensity/hadamard_code_mses', name, '-append')
% saveas(gcf, "Simulations-soft_thresholding/"+name, 'png');

function distanceImage = generateBaseImageGradientCircle(centerX,centerY,radius, dim1, dim2)
[columnsInImage, rowsInImage] = meshgrid(1:dim2, 1:dim1);

circlePixels = (rowsInImage - centerY).^2 + (columnsInImage - centerX).^2 <= radius.^2;
% Euclidean Distance Transform
distanceImage = bwdist(~circlePixels);
% Normalization
distanceImage = distanceImage / max(distanceImage(:));
distanceImage = double(distanceImage);
end