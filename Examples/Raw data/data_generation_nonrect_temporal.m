% Creates .mat file for personal tests on hadamard code, using generated
% Poisson optical sectioning microscopy patterns in complement interleaved
% fashion.
%

% tStart = tic;
seed = 0;

sparsityAdjustment = 56;

% image dimensions
npix_x = 600;
npix_y = 600;

nSeries = 1;

numPatterns1 = 16;
numPatterns2 = 16;

applyClosedLoopReconstruction(nSeries, numPatterns1, numPatterns1, sparsityAdjustment, npix_x, npix_y);

% tEnd = toc(tStart);

function finalImageSeries =  applyClosedLoopReconstruction(nSeries, numPatterns1, numPatterns2, sparsityAdjustment, npix_x, npix_y)
    finalImageSeries = zeros(npix_x, npix_y, nSeries);
    
    for seriesI = 1:nSeries
        %% PART 1 ------------> Run algorithm with fewer patterns

        % GENERATE ILLUMINATION PATTERNS, COMPLEMENT INTERLEAVED
        % Do so for each step in time series.

        % Define number of orthogonal locations.
        % Usually as n = m - 1 with m a multiple of 4,
        % in any case, m>n patterns are generated, with m as low as possible.
        % The offset defines how orthogonal locations distribute in space; they
        % should be equidistant from each other to minimize scattering crosstalk
        nlocations_and_offset = [(numPatterns1 - 1) 1]; % [n offset]

        npts = 120 - round(sparsityAdjustment / 2);
        sq_size = 56 - sparsityAdjustment;
        if (sq_size < 1 && sparsityAdjustment < 74)
            sq_size = 36 - round(sparsityAdjustment / 2) + 1;
        elseif (sq_size < 1)
            sq_size = 44 - round(sparsityAdjustment / 2);
        end
        spacing = sq_size * 2;

        % evenly spaced poisson disc patterns
        [poisson_patterns, without_complement] = generateEvenlySpacedPoissonDiscCodesMatrix(nlocations_and_offset, spacing, npts, npix_x, npix_y, sq_size);

        figure(22);
        moviesc(without_complement);
        title 'illumination patterns (evenly spaced)';

        % GENERATE BASE IMAGE (which should have same # pixels as in pattern)
        p = without_complement.data(:, :, 1) > 0;
        numPixels = sum(p, "all");
        [dim1, dim2, ~] = size(poisson_patterns);
        
        sparsity = ((npix_x * npix_y) - numPixels) / (npix_x * npix_y);
        disp(sparsity)

        distanceImage = generateNBinaryCircles(60, round(numPixels / nSeries), npix_x, npix_y); %3 --> nSeries

        figure(22);
        moviesc(vm(distanceImage));
        title 'image without noise';
        

        %%
        reconstructedImage = reconstructImage(distanceImage, poisson_patterns, without_complement);
        normalized_reconstruction = (reconstructedImage - min(reconstructedImage(:))) / (max(reconstructedImage(:) - min(reconstructedImage(:))));
        
        figure(23);
        moviesc(vm(normalized_reconstruction));
        title("Normalized Image Reconstructed using r patterns");

        %% PART 2  ------------> Determine ROI using reconstructed image, kmeans

        [clusterIndexes, clusterCenters] = kmeans(reconstructedImage(:), 2,...
          'distance', 'sqEuclidean', ...
          'Replicates', 2);
        clustering = reshape(clusterIndexes, dim1, dim2);
        
%         figure(300);
%         moviesc(vm(clustering));
%         title('Kmeans Clustering');

        disp(clusterCenters)
        [~, indexOfMaxValue] = max(clusterCenters);
        obj = clustering == indexOfMaxValue;

        props = regionprops(obj, 'FilledArea');
        props_arr = cell2mat(struct2cell(props)');

        % comparison # must be configured depending on image set
        num_clusters = sum(props_arr > .000000000000001);
        obj = bwareafilt(obj, num_clusters);

        [rows, columns, roi_vec] = find(obj .* normalized_reconstruction);

        %% PART 3 ------------> Run algorithm with more patterns using ROI
        roi = zeros(dim1, dim2);
        for i = 1:length(rows)
            roi(rows(i), columns(i)) = roi_vec(i);
        end
 
%         figure(202);
%         moviesc(vm(roi));

        % generate hadamard codes
        nlocations_and_offset = [(numPatterns2 - 1) 14]; % [n offset] % [(4 * index - 1) 1]
        % binning = [43 1]; % optional parameter binning to group projection elements

        [hadamard_patterns, without_complement] = generateEvenlySpacedPoissonDiscCodesVector(nlocations_and_offset, spacing, npts, rows, columns, sq_size);
        reconstructedImage2 = reconstructImageVectorized(roi_vec, hadamard_patterns, without_complement);

        %% PART 4 ------------> Reconstruct to original resolution of image

        [rdim1, rdim2, rdim3] = size(reconstructedImage);
        final = zeros(rdim1, rdim2, rdim3);
        for i = 1:length(rows)
            final(rows(i), columns(i)) = reconstructedImage2(i);
        end

        figure;
        moviesc(vm(final));
        title("final before tvd");
        
        finalImageSeries(:, :, seriesI) = reshape(im2double(mat2gray(applyTotalVariationDenoising(final))), dim1, dim2);
        
%         disp("MSE:");
%         disp(immse(distanceImage, finalImageSeries(:, :, seriesI)));
    end
    
    sparsity = ((npix_x * npix_y) - numPixels) / (npix_x * npix_y);
    disp(sparsity)
    
    figure(305);
    moviesc(vm(finalImageSeries));
    title("FINAL IMAGE SERIES WITH TVD");
end

function reconstructedImage = reconstructImage(distanceImage, patterns, without_complement)
    [dim1, dim2, dim3] = size(patterns);
    % ADD NOISE TO EACH FRAME, APPLY POISSON PATTERN
    [imgdim1, imgdim2, imgdim3] = size(distanceImage);
    reconstructedImage = zeros(imgdim1, imgdim2, imgdim3);
    for k = 1:imgdim3
    image_patterns = zeros(dim1, dim2, dim3);
    for i = 1:dim3
        pattern = patterns.data(:,:,i);
        new_img = distanceImage(:, :, k) .* pattern;
        new_img_with_gauss = imnoise(new_img,'gaussian');
        new_img_with_poiss = imnoise(new_img, 'poisson');
        image_patterns(:, :, i) = (new_img_with_gauss + new_img_with_poiss) / 2;
%         image_patterns(:, :, i) = new_img;
    end

    figure(4);
    new_image_patterns = vm(image_patterns);
    moviesc(new_image_patterns);


    % RECONSTRUCT IMAGE
    rootdir = fullfile('Examples/Raw data');
    res_fpath = fullfile(rootdir,'results.mat');

        %%
        ncomps = 1;
        %rng(0,'twister') % ensures multiple runs of eigs() will return repeatable results
        
        [sv, uhad, uref, uu] = dyn_had(new_image_patterns(:,:), without_complement(:,:), ncomps);
        %%
        save(res_fpath,...
            'sv', ...
            'uhad', ...
            'uref', ...
            'uu', ...
            'ncomps')
    %% 
    % extract time-integrated images  
    h50img = mean(vm(uhad*sv',new_image_patterns.imsz)); % hadamard image

    reconstructedImage(:, :, k) = h50img;
    end

    % mse_27_6_200_diff_patterns = immse(distanceImage, reconstructedImage);
    % name = 'mse_27_6_200_diff_patterns';
    % save('Simulations-mult_frames/hadamard_code_mses', name, '-append')
    % saveas(gcf, "Simulations-mult_frames/"+name, 'png');
end

function reconstructedImage = reconstructImageVectorized(distanceImage, hadamard_patterns, without_complement)
    [dim1, dim2] = size(hadamard_patterns);
    % ADD NOISE TO EACH FRAME, APPLY HADAMARD PATTERN
    [imgdim1, imgdim2, imgdim3] = size(distanceImage);
    reconstructedImage = zeros(imgdim1, imgdim2, imgdim3);
    for k = 1:imgdim3
    image_patterns = zeros(dim1, dim2);
    for i = 1:dim2
        pattern = hadamard_patterns.data(:, i);
        new_img = distanceImage(:, :, k) .* pattern;
        new_img_with_gauss = imnoise(new_img,'gaussian');
        new_img_with_poiss = imnoise(new_img, 'poisson');
        image_patterns(:, i) = (new_img_with_gauss + new_img_with_poiss) / 2;
%         image_patterns(:, i) = new_img;
    end

%     figure(4);
    new_image_patterns = vm(image_patterns);
%     moviesc(new_image_patterns);

    % RECONSTRUCT IMAGE
    rootdir = fullfile('Examples/Raw data');
    res_fpath = fullfile(rootdir,'results.mat');

    %%
    ncomps = 1;
    %rng(0,'twister') % ensures multiple runs of eigs() will return repeatable results

    [sv, uhad, uref, uu] = dyn_had_vec(new_image_patterns.data(:,:), without_complement.data(:,:), ncomps);
    %%
    save(res_fpath,...
        'sv', ...
        'uhad', ...
        'uref', ...
        'uu', ...
        'ncomps')
    %% 
    h50img = mean(uhad*sv', 2); % hadamard image
    h50img = vm(h50img);

    reconstructedImage(:, :, k) = h50img;
    end

    % mse_27_6_200_diff_patterns = immse(distanceImage, reconstructedImage);
    % name = 'mse_27_6_200_diff_patterns';
    % save('Simulations-mult_frames/hadamard_code_mses', name, '-append')
    % saveas(gcf, "Simulations-mult_frames/"+name, 'png');
end

function distanceImage = generateNBinaryCircles(n, totalPixels, dim1, dim2)
    [columnsInImage, rowsInImage] = meshgrid(1:dim2, 1:dim1);
    
    radius = round(sqrt(totalPixels / (n * pi)));
    
    centerX = randi(dim2, 1, n);
    centerY = randi(dim1, 1, n);
    circlePixels = (rowsInImage - centerY(1)).^2 + (columnsInImage - centerX(1)).^2 <= radius.^2;
    for circleI = 2:n
        tempPixels = (rowsInImage - centerY(circleI)).^2 + (columnsInImage - centerX(circleI)).^2 <= radius.^2;
        circlePixels = circlePixels | tempPixels;
    end
    
    % Normalization
    distanceImage = circlePixels / max(circlePixels(:));
    distanceImage = double(distanceImage)/10;

end

function distanceImage = generateNGradientCircles(n, totalPixels, dim1, dim2)
    [columnsInImage, rowsInImage] = meshgrid(1:dim2, 1:dim1);
    
    radius = round(sqrt(totalPixels / (n * pi)));
    
    centerX = randi(dim2, 1, n);
    centerY = randi(dim1, 1, n);
    circlePixels = (rowsInImage - centerY(1)).^2 + (columnsInImage - centerX(1)).^2 <= radius.^2;
    for circleI = 2:n
        tempPixels = (rowsInImage - centerY(circleI)).^2 + (columnsInImage - centerX(circleI)).^2 <= radius.^2;
        circlePixels = circlePixels | tempPixels;
    end

    % Euclidean Distance Transform
    distanceImage = bwdist(~circlePixels);
    % Normalization
    distanceImage = distanceImage / max(distanceImage(:));
    distanceImage = double(distanceImage)/10;

end

% function distanceImageSeries = generateTimeSeriesOfNBinaryCircles(nSeries, totalCircles, totalPixels, dim1, dim2)
%     numCirclesPerStep = round(totalCircles / nSeries);
%     totalPixelsPerStep = round(totalPixels / nSeries);
%     
%     distanceImageSeries = zeros(dim2, dim1, nSeries);
%     for seriesI = 1:nSeries
%         distanceImageSeries(:, :, seriesI) = generateNBinaryCircles(numCirclesPerStep, totalPixelsPerStep, dim1, dim2);
%     end
% end
% 
% function distanceImageSeries = generateTimeSeriesOfNGradientCircles(nSeries, totalCircles, totalPixels, dim1, dim2)
%     numCirclesPerStep = round(totalCircles / nSeries);
%     totalPixelsPerStep = round(totalPixels / nSeries);
%     
%     distanceImageSeries = zeros(dim2, dim1, nSeries);
%     for seriesI = 1:nSeries
%         distanceImageSeries(:, :, seriesI) = generateNGradientCircles(numCirclesPerStep, totalPixelsPerStep, dim1, dim2);
%     end
% end

function [disc_patterns, without_complement] = generateEvenlySpacedPoissonDiscCodesMatrix(nlocations_and_offset, pt_spacing, npts, npix_x, npix_y, sq_size)
num_frames = nlocations_and_offset(1) + 1;
patterns_logical = false(npix_x, npix_y, num_frames);

% numberPts = round((npix_x * npix_y) / ((pt_spacing / 3) ^ 2)) * num_frames;
% [all_pts] = round(poissonDisc([npix_x, npix_y], pt_spacing / num_frames, numberPts, 0));

gridNumX = ceil(npix_x / (pt_spacing / 2));
gridNumY = ceil(npix_y / (pt_spacing / 2));
[all_pts] = round(poissonDisc([gridNumX, gridNumY], 0.1));
all_pts = unique(all_pts * pt_spacing, 'rows');

all_pts_rep = repelem(all_pts, 5, 1);

% [all_pts] = round(poissonDisc([gridNumX, gridNumY, num_frames], 1));
% all_pts(:, 1, :) = all_pts(:, 1, :) * pt_spacing;
% all_pts(:, 2, :) = all_pts(:, 2, :) * pt_spacing;
% all_pts = unique(all_pts(:, 1:2, :), 'rows');

num_all_pts = length(all_pts_rep);

for j=1:num_frames
    % set random seed
%     rng(j * 10);
%     
    sq_size_multiplier = 1 + 1.25*rand();
    disp(sq_size_multiplier);
    sq_size = sq_size * sq_size_multiplier;  
    
    indices = randsample(1:length(all_pts), floor(num_all_pts / num_frames), false); 
    pts = all_pts(indices, :);
    
    toRemove = pts;
    for r=1:length(toRemove)
        idx = find(all_pts_rep(:, 1) == toRemove(r, 1) & ...
            all_pts_rep(:, 2) == toRemove(r, 2), 1);
        all_pts_rep(idx, :) = [];
    end
    
    if j ~= num_frames - 1
        all_pts = unique(all_pts_rep, 'rows');
    end

    pts(pts(:,1) < sq_size+1,:)=[];
    pts(pts(:,2) < sq_size+1,:)=[];
    pts(pts(:,2) > 600 - sq_size,:)=[];
    pts(pts(:,1) > 600 - sq_size,:)=[];
    
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

figure(21);
moviesc(vm(sum(without_complement.data, 3)));
title 'sum of r illumination patterns (evenly spaced)';
end

function [disc_patterns, without_complement] = generateEvenlySpacedPoissonDiscCodesVector(nlocations_and_offset, pt_spacing, npts, index_rows, index_cols, sq_size)
    num_frames = nlocations_and_offset(1) + 1;
    
    roi_rows = ceil((max(index_rows) - min(index_rows)) * 1.15);
    roi_cols = ceil((max(index_cols) - min(index_cols)) * 1.15);
    
    patterns_logical = false(roi_rows, roi_cols, num_frames);
    
% %     numberPts = round(length(index_rows) / ((pt_spacing / 3))) * num_frames;
%     numberPts = round(length(index_rows) / ((pt_spacing / 3) ^ 2)) * num_frames;
%     [all_pts] = round(poissonDisc([roi_rows, roi_cols], pt_spacing / num_frames, numberPts, 0));

    gridNumX = ceil(roi_rows / (pt_spacing / 2));
    gridNumY = ceil(roi_cols / (pt_spacing / 2));
    [all_pts] = round(poissonDisc([gridNumX, gridNumY], 0.1));
    [all_pts] = unique(all_pts * pt_spacing, 'rows');
    
    all_pts_rep = repelem(all_pts, 5, 1);
    
    num_all_pts = length(all_pts_rep);

    for j=1:num_frames
%         % set random seed
%         rng(j * 10);

        sq_size_multiplier = 1 + 1.25*rand();
        disp(sq_size_multiplier);
        sq_size = floor(sq_size * sq_size_multiplier);
    
        indices = randsample(1:length(all_pts), floor(num_all_pts / num_frames), false); 
        pts = all_pts(indices, :);

        toRemove = pts;
        for r=1:length(toRemove)
            idx = find(all_pts_rep(:, 1) == toRemove(r, 1) & ...
                all_pts_rep(:, 2) == toRemove(r, 2), 1);
            all_pts_rep(idx, :) = [];
        end

        if j ~= num_frames - 1
            all_pts = unique(all_pts_rep, 'rows');
        end

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
%     figure(1)
%     moviesc(patterns)
%     title 'illumination patterns, second round'
% 
%     figure(2)
%     moviesc(without_comp)
%     title 'illumination patterns without interleaved complements, second round'
%     
    figure(3)
    moviesc(vm(sum(without_comp.data, 3)))
    title 'sum of illumination patterns (evenly spaced)'
    
    
    % Reshape 3d matrices to 2d matrices where each column is an image 
    % frame vector
    hdim3 = size(patterns, 3);
    wdim3 = size(without_comp, 3);
    vec_size = length(index_rows);
    
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
    img = SB_ATV(reconstructedImg, 0.25);
end