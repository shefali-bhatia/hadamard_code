% Creates .mat file for personal tests on hadamard code, using generated
% Hadamard optical sectioning microscopy patterns in complement interleaved
% fashion.
%

% CLOSED LOOP PATTERN
% PART 1 ------------> Run algorithm with fewer patterns

% GENERATE HADAMARD CODES, COMPLEMENT INTERLEAVED

% Define number of orthogonal locations.
% Usually as n = m - 1 with m a multiple of 4,
% in any case, m>n patterns are generated, with m as low as possible.
% The offset defines how orthogonal locations distribute in space; they
% should be equidistant from each other to minimize scattering crosstalk
nlocations_and_offset = [19 5]; % [n offset]

% optional parameter binning to group projection elements
binning = [150 150];
spacing = 75;
npts = 120;
sq_size = 55;

% image dimensions
npix_x = 1200;
npix_y = 1200;

[hadamard_patterns, without_complement] = generatePoissonDiscCodesMatrix(nlocations_and_offset, spacing, npts, npix_x, npix_y, sq_size);

% disp(min(hadamard_patterns.data(:)))
% disp(max(hadamard_patterns.data(:)))
% Display patterns as imagej stack style. note there are 2*m frames
figure(1)
moviesc(hadamard_patterns)
title 'illumination patterns'

figure(2)
moviesc(without_complement)
title 'illumination patterns without interleaved complements'

% GENERATE BASE IMAGE
[dim1, dim2, dim3] = size(hadamard_patterns);

centerX1 = 600;
centerY1 = 600;
centerX2 = 800;
centerY2 = 600;
radius = 400;

centerX = npix_x/2;
centerY = npix_y/2;

% distanceImage = generate2BaseImageGradientCircles(centerX1,centerY1,centerX2,centerY2,radius, dim1, dim2);
distanceImage = generateBaseImageGradientCircle(centerX,centerY,radius, dim1, dim2);

figure(3);
d = vm(distanceImage);
moviesc(d);

% GENERATE "MOVING" BASE IMAGES
% [dim1, dim2, dim3] = size(hadamard_patterns);
% 
% timestep1 = generateBaseImageGradientCircle(20,22,10,dim1,dim2);
% timestep2 = generateBaseImageGradientCircle(40,22,10,dim1,dim2);
% timestep3 = generateBaseImageGradientCircle(60,22,10,dim1,dim2);
% 
% distanceImage = cat(3, timestep1, timestep2, timestep3);
% 
% figure(3);
% d = vm(distanceImage);
% moviesc(d);

reps = 1;

reconstructedImage = reconstructImage(distanceImage, hadamard_patterns, without_complement, reps);

figure(200);
moviesc(vm(reconstructedImage));
title("Image Reconstructed using r patterns");

% PART 2  ------------> Determine ROI using reconstructed image, kmeans

mask = reconstructedImage & (reconstructedImage > 1e-5);
figure(305);
moviesc(vm(mask));

[clusterIndexes, clusterCenters] = kmeans(reconstructedImage(:), 2,...
  'distance', 'sqEuclidean', ...
  'Replicates', 2);
clustering = reshape(clusterIndexes, dim1, dim2);
figure(300);
moviesc(vm(clustering));
title('Kmeans Clustering');

% FIXME --> change to use greythresh, find(distanceImage > T)
%   --> nvm that did not work
disp(clusterCenters)
[~, indexOfMaxValue] = max(clusterCenters);
obj = clustering == indexOfMaxValue;
% obj = imbinarize(obj .* reconstructedImage, graythresh(obj .* reconstructedImage));


props = regionprops(obj, 'FilledArea');
props_arr = cell2mat(struct2cell(props)');
max_prop = max(props_arr);
i = 1;
for j = 1:length(props)
    if props_arr(i)/max_prop < .5
        props_arr(i) = [];
    else
        i = i + 1;
    end
end

num_clusters = length(props_arr);
obj = bwareafilt(obj, num_clusters);


figure(301);
moviesc(vm(obj));
title('Max Cluster');

% % indices = find(obj);
% % [r_ind, c_ind] = ind2sub([npix_y, npix_x], indices);
% % disp([r_ind, c_ind]);

% extract image roi
% % r_range = [min(r_ind), max(r_ind)];
% % disp(r_range);
% % c_range = [min(c_ind), max(c_ind)];
% % disp(c_range);

[rows, columns, roi_vec] = find(obj .* distanceImage);


% PART 3 ------------> Run algorithm with more patterns using ROI
% % roi = distanceImage(r_range(1):r_range(2), c_range(1):c_range(2));
roi = zeros(dim1, dim2);
for i = 1:length(rows)
    roi(rows(i), columns(i)) = roi_vec(i);
end
figure(202);
moviesc(vm(roi));

% generate hadamard codes
nlocations_and_offset = [63 14]; % [n offset]
binning = [43 1]; % optional parameter binning to group projection elements

% % rdim = r_range(2)-r_range(1)+1;
% % cdim = c_range(2)-c_range(1)+1;

[hadamard_patterns, without_complement] = generatePoissonDiscCodesVector(nlocations_and_offset, spacing, npts, rows, columns, sq_size);

% Display patterns as imagej stack style. note there are 2*m frames
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


% figure(1)
% moviesc(vm(hadamard_patterns_matrix))
% title 'illumination patterns, second round'
% 
% figure(2)
% moviesc(vm(without_complement_matrix))
% title 'illumination patterns without interleaved complements, second round'

reps = 1;

reconstructedImage2 = reconstructImageVectorized(roi_vec, hadamard_patterns, without_complement, reps);

% % use indices from find operation to transform vector back into image dim
% 
% 
% 
% figure(203);
% moviesc(vm(reconstructedImage2));
% title("Second round of reconstruction");

% PART 4 ------------> Reconstruct to original resolution of image

[rdim1, rdim2, rdim3] = size(reconstructedImage);
final = zeros(rdim1, rdim2, rdim3);
% final(r_range(1):r_range(2), c_range(1):c_range(2)) = reconstructedImage2(:,:);
% % final = final .* obj;
% % final(final<0) = 0;
for i = 1:length(rows)
    final(rows(i), columns(i)) = reconstructedImage2(i);
end


figure(204);
moviesc(vm(final));
title("FINAL IMAGE");

final_with_tvd = reshape(applyTotalVariationDenoising(final), dim1, dim2);
figure(205);
moviesc(vm(final_with_tvd));
title("FINAL IMAGE WITH TVD");


function reconstructedImage = reconstructImage(distanceImage, hadamard_patterns, without_complement, reps)
    [dim1, dim2, dim3] = size(hadamard_patterns);
    % ADD NOISE TO EACH FRAME, APPLY HADAMARD PATTERN
    [imgdim1, imgdim2, imgdim3] = size(distanceImage);
    reconstructedImage = zeros(imgdim1, imgdim2, imgdim3);
    for k = 1:imgdim3
    image_patterns = zeros(dim1, dim2, dim3*reps);
    disp(k)
    for i = 1:(dim3*reps)
        j = mod(i, dim3);
        if j == 0
            j = dim3;
        end
        pattern = hadamard_patterns.data(:,:,j);
        new_img = distanceImage(:, :, k) .* pattern;
%         new_img_with_gauss = imnoise(new_img,'gaussian');
%         new_img_with_poiss = imnoise(new_img, 'poisson');
%         image_patterns(:, :, i) = new_img_with_gauss + new_img_with_poiss;
        image_patterns(:, :, i) = new_img;
    end

    figure(4);
    new_image_patterns = vm(image_patterns);
    moviesc(new_image_patterns);



    % RECONSTRUCT IMAGE
    rootdir = fullfile('Examples/Raw data');
    %% Code parameters
    il = 20*2; % interleaved length containing 2 series of length 20
    hl = il/2; % hadamard length, equal to half the interleaved length
    hadtraces = hadamard_bincode_nopermutation(hl-1)'*2-1;

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
        [sv, uhad, uref, uu] = dyn_had(new_image_patterns(:,:), without_complement(:,:), ncomps);
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

    reconstructedImage(:, :, k) = h50img;

    % figure(100+k);
    % moviesc(vm(h50img));
    end

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
        [sv, uhad, uref, uu] = dyn_had_vec(new_image_patterns.data(:,:), without_complement.data(:,:), ncomps);
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

function distanceImage = generateBaseImageGradientCircle(centerX,centerY,radius, dim1, dim2)
    [columnsInImage, rowsInImage] = meshgrid(1:dim2, 1:dim1);

    circlePixels = (rowsInImage - centerY).^2 + (columnsInImage - centerX).^2 <= radius.^2;
    % Euclidean Distance Transform
    distanceImage = bwdist(~circlePixels);
    % Normalization
    distanceImage = distanceImage / max(distanceImage(:));
    distanceImage = double(distanceImage)/50;
end

function distanceImage = generate2BaseImageGradientCircles(centerX1,centerY1,centerX2, centerY2, radius, dim1, dim2)
    [columnsInImage, rowsInImage] = meshgrid(1:dim2, 1:dim1);

    circle1Pixels = (rowsInImage - centerY1).^2 + (columnsInImage - centerX1).^2 <= radius.^2;
    circle2Pixels = (rowsInImage - centerY2).^2 + (columnsInImage - centerX2).^2 <= radius.^2;
    circlePixels = circle1Pixels | circle2Pixels;
    % Euclidean Distance Transform
    distanceImage = bwdist(~circlePixels);
    % Normalization
    distanceImage = distanceImage / max(distanceImage(:));
    distanceImage = double(distanceImage)/100;
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

function [hadamard_patterns, without_complement] = generateHadamardCodesMatrix(nlocations_and_offset, binning, npix_x, npix_y)% Generate a matrix of size [1024 768]. randomization is always the same.
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
hadamard_patterns = hadamard_patterns(1:npix_y,1:npix_x,:);
without_complement = without_complement(1:npix_y,1:npix_x,:);
end

function [rand_patterns, without_complement] = generateRandomCodesMatrix(nlocations_and_offset, binning, npix_x, npix_y)
patterns_logical_bins = logical(randi([0, 1], npix_x / binning(1), npix_y  / binning(2), nlocations_and_offset(1) + 1));

[x, y, z] = size(patterns_logical_bins);

patterns_logical = false(npix_x, npix_y, z);
for p=1:z
    patterns_logical(:, :, p) = imresize(patterns_logical_bins(:, :, p), [x y] .* binning, 'nearest');
end

% Interleave with complement in 3rd dimension, store as vectorized movie
toint = {... % to interleave matrices in 3d, place them in a cell to then use the permute-cell2mat-permute-reshape method
     patterns_logical ,...
    ~patterns_logical};
rand_patterns = vm(reshape(permute(cell2mat(permute(toint,[1 4 3 2])),[1 2 4 3]),size(toint{1},1),size(toint{1},2),[]));
without_complement = vm(patterns_logical);
end

function [rand_patterns, without_complement] = generateRandomCodesVector(nlocations_and_offset, binning, index_rows, index_cols)
    patterns_logical_bins = zeros(length(index_rows) / binning(1), 1, nlocations_and_offset(1) + 1);
    for i=1:nlocations_and_offset(1) + 1
        patterns_logical_bins(:,:,i) = logical(randi([0, 1], length(index_rows) / binning(1), length(index_cols)));
    end
    
    [x, y, z] = size(patterns_logical_bins);

    patterns_logical = false(length(index_rows) , length(index_cols), z);
    for p=1:z
        patterns_logical(:, :, p) = imresize(patterns_logical_bins(:, :, p), [x y] .* binning, 'nearest');
    end

    % Interleave with complement in 3rd dimension, store as vectorized movie
    toint = {... % to interleave matrices in 3d, place them in a cell to then use the permute-cell2mat-permute-reshape method
         patterns_logical ,...
        ~patterns_logical};
    patterns = vm(reshape(permute(cell2mat(permute(toint,[1 4 3 2])),[1 2 4 3]),size(toint{1},1),size(toint{1},2),[]));
    without_comp = vm(patterns_logical);

    hdim3 = size(patterns, 3);
    wdim3 = size(without_comp, 3);
    vec_size = length(index_rows);
    
    rand_patterns = zeros(vec_size, hdim3);
    without_complement = zeros(vec_size, wdim3);
    
    for d = 1:hdim3
        for i = 1:vec_size
            r = index_rows(i);
            c = index_cols(i);
            rand_patterns(i, d) = patterns(r, c, d);
        end
    end
    
    for d = 1:wdim3
        for i = 1:vec_size
            r = index_rows(i);
            c = index_cols(i);
            without_complement(i, d) = without_comp(r, c, d);
        end
    end
    rand_patterns = vm(rand_patterns);
    without_complement = vm(without_complement);
end

function [disc_patterns, without_complement] = generatePoissonDiscCodesMatrix(nlocations_and_offset, pt_spacing, npts, npix_x, npix_y, sq_size)
num_frames = nlocations_and_offset(1) + 1;
patterns_logical = false(npix_x, npix_y, num_frames);

for j=1:num_frames+1
    [pts] = round(poissonDisc([npix_x, npix_y], pt_spacing, npts, 0));

    pts(pts(:,1) < sq_size+1,:)=[];
    pts(pts(:,2) < sq_size+1,:)=[];
    pts(pts(:,2) > 1200 - sq_size,:)=[];
    pts(pts(:,1) > 1200 - sq_size,:)=[];
    
    for i = 1:size(pts, 1)
        patterns_logical(pts(i, 1)-sq_size:pts(i, 1)+sq_size, pts(i, 2)-sq_size:pts(i, 2)+sq_size, j) = 1;
    end
end

% Interleave with complement in 3rd dimension, store as vectorized movie
toint = {... % to interleave matrices in 3d, place them in a cell to then use the permute-cell2mat-permute-reshape method
     patterns_logical ,...
    ~patterns_logical};
disc_patterns = vm(reshape(permute(cell2mat(permute(toint,[1 4 3 2])),[1 2 4 3]),size(toint{1},1),size(toint{1},2),[]));
without_complement = vm(patterns_logical);
end

function [disc_patterns, without_complement] = generatePoissonDiscCodesVector(nlocations_and_offset, pt_spacing, npts, index_rows, index_cols, sq_size)
    num_frames = nlocations_and_offset(1) + 1;
    
    roi_rows = ceil((max(index_rows) - min(index_rows)) * 1.15);
    roi_cols = ceil((max(index_cols) - min(index_cols)) * 1.15);
    
    patterns_logical = false(roi_rows, roi_cols, num_frames);

    for j=1:num_frames+1
        [pts] = round(poissonDisc([roi_rows, roi_cols], pt_spacing, npts, 0));

        pts(pts(:,1) < sq_size+1,:)=[];
        pts(pts(:,2) < sq_size+1,:)=[];
        pts(pts(:,2) > roi_cols - sq_size,:)=[];
        pts(pts(:,1) > roi_cols - sq_size,:)=[];

        for i = 1:size(pts, 1)
            patterns_logical(pts(i, 1)-sq_size:pts(i, 1)+sq_size, pts(i, 2)-sq_size:pts(i, 2)+sq_size, j) = 1;
        end
        
        % code for poisson disc circle instead of square
%         circlePts = (patterns_logical(:, :, j) - pts(:, 1)).^2 + (patterns_logical(:, :, j) - pts(:, 2)).^2 <= sq_size.^2;
%         distanceImage = bwdist(~circlePts);
%         patterns_logical(pts(i, 1), pts(i, 2), j) = 1;
        
        
%         figure();
%         imagesc(patterns_logical(:,:,j));
        
    end

    % Interleave with complement in 3rd dimension, store as vectorized movie
    toint = {... % to interleave matrices in 3d, place them in a cell to then use the permute-cell2mat-permute-reshape method
         patterns_logical ,...
        ~patterns_logical};
    patterns = vm(reshape(permute(cell2mat(permute(toint,[1 4 3 2])),[1 2 4 3]),size(toint{1},1),size(toint{1},2),[]));
    without_comp = vm(patterns_logical);
    
    
    % Display patterns
    figure(1)
    moviesc(patterns)
    title 'illumination patterns, second round'

    figure(2)
    moviesc(without_comp)
    title 'illumination patterns without interleaved complements, second round'
    
    
    % Reshape 3d matrices to 2d matrices where each column is an image 
    % frame vector
    hdim3 = size(patterns, 3);
    wdim3 = size(without_comp, 3);
    vec_size = length(index_rows);
    
%     disc_patterns = zeros(vec_size, hdim3);
%     without_complement = zeros(vec_size, wdim3);
    
    r = index_rows(:) - min(index_rows) + 1;
    c = index_cols(:) - min(index_cols) + 1;
    
    R = repmat(r, [hdim3, 1]);
    C = repmat(c, [hdim3, 1]);
    
    K = kron((1:hdim3)', ones(numel(r), 1));
    ind = sub2ind(size(patterns),R,C,K);
    
    disc_patterns = reshape(patterns.data(ind), vec_size, hdim3);
    
    
    R2 = repmat(r, [wdim3, 1]);
    C2 = repmat(c, [wdim3, 1]);
    
    K2 = kron((1:wdim3)', ones(numel(r), 1));
    ind2 = sub2ind(size(patterns),R2, C2, K2);
    
    without_complement = reshape(without_comp.data(ind2), vec_size, wdim3);
    
    disc_patterns = vm(disc_patterns);
    without_complement = vm(without_complement);
end

function img = applyTotalVariationDenoising(reconstructedImg)
    img = SB_ATV(reconstructedImg, 0.5);
end