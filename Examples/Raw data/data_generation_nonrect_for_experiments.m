% Creates .mat file for personal tests on hadamard code, using generated
% Hadamard optical sectioning microscopy patterns in complement interleaved
% fashion.
%

% Dimensions of images taken during experiment
dim1 = 600;
dim2 = 600;

% LOAD IMAGES THAT WERE TAKEN DURING EXPERIMENT
filename = fullfile('SI/3.9umRedbeads_region1_1/Pos0', '*.tif');
directory = dir(filename);

% last image in the directory is not applicable
image_patterns = zeros(dim1, dim2, length(directory)-1);
for k=1:length(directory)-1
    file = strcat('SI/3.9umRedbeads_region1_1/Pos0/', directory(k).name);
    image = mat2gray(imread(file));
    image_patterns(:, :, k) = im2double(image);
end

figure(1);
moviesc(vm(image_patterns));
title('images taken using interpolated hadamard patterns');

% LOAD INTERPOLATED HADAMARD PATTERNS
filename=fullfile('SI/homogeneousSample_1/Pos0', '*.tif');
directory = dir(filename);

% last image in the directory is not applicable
interpolated_patterns = zeros(dim1, dim2, length(directory)-1);
bin_interpolated_patterns = zeros(dim1, dim2, length(directory)-1);
for k=1:length(directory)-1
    file = strcat('SI/homogeneousSample_1/Pos0/', directory(k).name);
    interpolated_patterns(:, :, k) = imread(file);
    bin_interpolated_patterns(:, :, k) = imbinarize(imread(file));
end

figure(2);
moviesc(vm(interpolated_patterns));
title('interpolated hadamard patterns');

interp_without_complement = bin_interpolated_patterns(:,:,1:2:end);

figure(3);
moviesc(vm(interp_without_complement));
title('interpolated hadamard patterns without complements');

figure(4);
moviesc(vm(bin_interpolated_patterns));
title('interpolated hadamard patterns binary');

% LOAD NON-INTERPOLATED HADAMARD PATTERNS
load('2560x1600_hadamard_patterns.mat');

had_patterns = hadamard_pattern_data(1:dim1, 1:dim2, :);
hadamard_pattern_without_comp = had_patterns(:,:,1:2:end);


figure(5);
moviesc(vm(hadamard_pattern_data));
title('generated hadamard patterns');

figure(6);
moviesc(vm(hadamard_pattern_without_comp));
title('generated hadamard patterns');


% CLOSED LOOP PATTERN
% PART 1 ------------> Run algorithm with fewer patterns
reconstructedImage = reconstructImage(image_patterns, interp_without_complement);


figure(200);
moviesc(vm(reconstructedImage));
title("Image Reconstructed using r patterns");

% PART 2  ------------> Determine ROI using reconstructed image, kmeans

% mask = reconstructedImage & (reconstructedImage > 1e-5);
% figure(305);
% moviesc(vm(mask));
% 
% [clusterIndexes, clusterCenters] = kmeans(reconstructedImage(:), 2,...
%   'distance', 'sqEuclidean', ...
%   'Replicates', 2);
% clustering = reshape(clusterIndexes, 44, 80);
% figure(300);
% moviesc(vm(clustering));
% title('Kmeans Clustering');
% 
% % FIXME --> change to use greythresh, find(distanceImage > T)
% %   --> nvm that did not work
% disp(clusterCenters)
% [~, indexOfMaxValue] = max(clusterCenters);
% obj = clustering == indexOfMaxValue;
% % obj = imbinarize(obj .* reconstructedImage, graythresh(obj .* reconstructedImage));
% 
% 
% props = regionprops(obj, 'FilledArea');
% props_arr = cell2mat(struct2cell(props)');
% max_prop = max(props_arr);
% i = 1;
% for j = 1:length(props)
%     if props_arr(i)/max_prop < .5
%         props_arr(i) = [];
%     else
%         i = i + 1;
%     end
% end
% 
% num_clusters = length(props_arr);
% obj = bwareafilt(obj, num_clusters);
% 
% 
% figure(301);
% moviesc(vm(obj));
% title('Max Cluster');
% 
% % % indices = find(obj);
% % % [r_ind, c_ind] = ind2sub([npix_y, npix_x], indices);
% % % disp([r_ind, c_ind]);
% 
% % extract image roi
% % % r_range = [min(r_ind), max(r_ind)];
% % % disp(r_range);
% % % c_range = [min(c_ind), max(c_ind)];
% % % disp(c_range);
% 
% [rows, columns, roi_vec] = find(obj .* distanceImage);
% 
% 
% % PART 3 ------------> Run algorithm with more patterns using ROI
% % % roi = distanceImage(r_range(1):r_range(2), c_range(1):c_range(2));
% roi = zeros(dim1, dim2);
% for i = 1:length(rows)
%     roi(rows(i), columns(i)) = roi_vec(i);
% end
% figure(202);
% moviesc(vm(roi));
% 
% % generate hadamard codes
% nlocations_and_offset = [63 14]; % [n offset]
% binning = [1 1]; % optional parameter binning to group projection elements
% 
% % % rdim = r_range(2)-r_range(1)+1;
% % % cdim = c_range(2)-c_range(1)+1;
% 
% [hadamard_patterns, without_complement] = generateHadamardCodesVector(nlocations_and_offset, binning, rows, columns);
% 
% % Display patterns as imagej stack style. note there are 2*m frames
% hadamard_patterns_matrix = zeros(dim1, dim2, size(hadamard_patterns, 2));
% without_complement_matrix = zeros(dim1, dim2, size(without_complement, 2));
% for d = 1:size(hadamard_patterns, 2)
%     for i = 1:length(rows)
%         hadamard_patterns_matrix(rows(i), columns(i), d) = hadamard_patterns.data(i, d);
%     end
% end
% 
% for d = 1:size(without_complement, 2)
%     for i = 1:length(rows)
%         without_complement_matrix(rows(i), columns(i), d) = without_complement.data(i, d);
%     end
% end
% 
% 
% figure(1)
% moviesc(vm(hadamard_patterns_matrix))
% title 'illumination patterns, second round'
% 
% figure(2)
% moviesc(vm(without_complement_matrix))
% title 'illumination patterns without interleaved complements, second round'
% 
% reps = 240;
% 
% reconstructedImage2 = reconstructImageVectorized(roi_vec, hadamard_patterns, without_complement, reps);
% 
% % % use indices from find operation to transform vector back into image dim
% % 
% % 
% % 
% % figure(203);
% % moviesc(vm(reconstructedImage2));
% % title("Second round of reconstruction");
% 
% % PART 4 ------------> Reconstruct to original resolution of image
% 
% [rdim1, rdim2, rdim3] = size(reconstructedImage);
% final = zeros(rdim1, rdim2, rdim3);
% % final(r_range(1):r_range(2), c_range(1):c_range(2)) = reconstructedImage2(:,:);
% % % final = final .* obj;
% % % final(final<0) = 0;
% for i = 1:length(rows)
%     final(rows(i), columns(i)) = reconstructedImage2(i);
% end
% 
% 
% figure(204);
% moviesc(vm(final));
% title("FINAL IMAGE");


function reconstructedImage = reconstructImage(image_patterns, without_complement)
    [imgdim1, imgdim2, ~] = size(image_patterns);
    reconstructedImage = zeros(imgdim1, imgdim2);

    new_image_patterns = vm(image_patterns);

    % RECONSTRUCT IMAGE
    rootdir = fullfile('Examples/Raw data');
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
    w50img = mean(vm(uref*sv',new_image_patterns.imsz)); % widefield image
    h50img = mean(vm(uhad*sv',new_image_patterns.imsz)); % hadamard image

    % reload reconstructed data from depth 60 um
    % res_fpath = fullfile(rootdir,'results.mat');
    % load(res_fpath,...
    %     'sv', ...
    %     'uhad', ...
    %     'uref', ...
    %     'uu', ...
    %     'ncomps')
    % % extract time-integrated images 
    % w60img = mean(vm(uref*sv',image_patterns.imsz));
    % h60img = mean(vm(uhad*sv',image_patterns.imsz));

    reconstructedImage(:, :) = h50img;

    % figure(100+k);
    % moviesc(vm(h50img));

    % mse_27_6_200_diff_patterns = immse(distanceImage, reconstructedImage);
    % name = 'mse_27_6_200_diff_patterns';
    % save('Simulations-mult_frames/hadamard_code_mses', name, '-append')
    % saveas(gcf, "Simulations-mult_frames/"+name, 'png');
end

function reconstructedImage = reconstructImageVectorized(distanceImage, hadamard_patterns, without_complement, reps)
    [dim1, dim2] = size(hadamard_patterns);
    % ADD NOISE TO EACH FRAME, APPLY HADAMARD PATTERN
    [imgdim1, imgdim2, imgdim3] = size(distanceImage);
    reconstructedImage = zeros(imgdim1, imgdim2, imgdim3);
    for k = 1:imgdim3
    image_patterns = zeros(dim1, dim2*reps);
    disp(k)
    for i = 1:(dim2*reps)
        j = mod(i, dim2);
        if j == 0
            j = dim2;
        end
        pattern = hadamard_patterns.data(:, j);
        new_img = distanceImage(:, :, k) .* pattern;
        new_img_with_gauss = imnoise(new_img,'gaussian');
        new_img_with_poiss = imnoise(new_img, 'poisson');
        image_patterns(:, i) = new_img_with_gauss + new_img_with_poiss;
    end

%     figure(4);
    new_image_patterns = vm(image_patterns);
%     moviesc(new_image_patterns);



    % RECONSTRUCT IMAGE
    rootdir = fullfile('Examples/Raw data');
    %% Code parameters
%     il = 20*2; % interleaved length containing 2 series of length 20
%     hl = il/2; % hadamard length, equal to half the interleaved length
%     hadtraces = hadamard_bincode_nopermutation(hl-1)'*2-1;

    %% Format calibration data
    % horiz_dim = 80;
    % vert_dim = 44;
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
        [sv, uhad, uref, uu] = dyn_had_vec(new_image_patterns.data(:,:), without_complement.data(:,:), ncomps, dim1, dim2);
        %%
        save(res_fpath,...
            'sv', ...
            'uhad', ...
            'uref', ...
            'uu', ...
            'ncomps')
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
%     w50img = mean(vm(uref*sv',new_image_patterns.imsz)); % widefield image
    h50img = mean(uhad*sv', 2); % hadamard image
    h50img = vm(h50img);

    % reload reconstructed data from depth 60 um
    % res_fpath = fullfile(rootdir,'results.mat');
    % load(res_fpath,...
    %     'sv', ...
    %     'uhad', ...
    %     'uref', ...
    %     'uu', ...
    %     'ncomps')
    % % extract time-integrated images 
    % w60img = mean(vm(uref*sv',image_patterns.imsz));
    % h60img = mean(vm(uhad*sv',image_patterns.imsz));

    reconstructedImage(:, :, k) = h50img;

    % figure(100+k);
    % moviesc(vm(h50img));
    end

    % mse_27_6_200_diff_patterns = immse(distanceImage, reconstructedImage);
    % name = 'mse_27_6_200_diff_patterns';
    % save('Simulations-mult_frames/hadamard_code_mses', name, '-append')
    % saveas(gcf, "Simulations-mult_frames/"+name, 'png');
end

function [hadamard_patterns, without_complement] = generateHadamardCodesVector(nlocations_and_offset, binning, index_rows, index_cols)
    % Generate a matrix of size [1024 768]. randomization is always the same.
    patterns_logical = alp_btd_to_logical(hadamard_patterns_scramble_nopermutation(nlocations_and_offset,binning));

    % Interleave with complement in 3rd dimension, store as vectorized movie
    toint = {... % to interleave matrices in 3d, place them in a cell to then use the permute-cell2mat-permute-reshape method
         patterns_logical ,...
        ~patterns_logical};
    hpatterns = vm(reshape(permute(cell2mat(permute(toint,[1 4 3 2])),[1 2 4 3]),size(toint{1},1),size(toint{1},2),[]));
    without_comp = vm(patterns_logical);

    % figure(5)
    % moviesc(hadamard_patterns)

    % Crop the matrix around ROI for easier visualization
    % hadamard_patterns = hadamard_patterns(1:npix_y,1:npix_x,:);
    % without_complement = without_complement(1:npix_y,1:npix_x,:);
    hdim3 = size(hpatterns, 3);
    wdim3 = size(without_comp, 3);
    vec_size = length(index_rows);
    
    hadamard_patterns = zeros(vec_size, hdim3);
    without_complement = zeros(vec_size, wdim3);
    
    for d = 1:hdim3
        for i = 1:vec_size
            r = index_rows(i);
            c = index_cols(i);
            hadamard_patterns(i, d) = hpatterns(r, c, d);
        end
    end
    
    for d = 1:wdim3
        for i = 1:vec_size
            r = index_rows(i);
            c = index_cols(i);
            without_complement(i, d) = without_comp(r, c, d);
        end
    end
    hadamard_patterns = vm(hadamard_patterns);
    without_complement = vm(without_complement);
end

% function [hadamard_patterns, without_complement] = generateHadamardCodesMatrix(nlocations_and_offset, binning, npix_x, npix_y)% Generate a matrix of size [1024 768]. randomization is always the same.
% patterns_logical = alp_btd_to_logical(hadamard_patterns_scramble_nopermutation(nlocations_and_offset,binning));
% 
% % Interleave with complement in 3rd dimension, store as vectorized movie
% toint = {... % to interleave matrices in 3d, place them in a cell to then use the permute-cell2mat-permute-reshape method
%      patterns_logical ,...
%     ~patterns_logical};
% hadamard_patterns = vm(reshape(permute(cell2mat(permute(toint,[1 4 3 2])),[1 2 4 3]),size(toint{1},1),size(toint{1},2),[]));
% without_complement = vm(patterns_logical);
% 
% % figure(5)
% % moviesc(hadamard_patterns)
% 
% % Crop the matrix around ROI for easier visualization
% hadamard_patterns = hadamard_patterns(1:npix_y,1:npix_x,:);
% without_complement = without_complement(1:npix_y,1:npix_x,:);
% end