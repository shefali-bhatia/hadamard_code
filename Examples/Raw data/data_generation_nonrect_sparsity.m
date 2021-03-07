% Creates .mat file for personal tests on hadamard code, using generated
% Poisson optical sectioning microscopy patterns in complement interleaved
% fashion.
%

%% Number of Patterns Loop
pArray = [4 8 16];
for numPatterns = 1:3
    p_index = pArray(numPatterns);
    legend_labels = [];
    %% Seed Number Loop
    for index = 1:5
        rng(index * 10);
        disp(rng);

        %% Simulation Code Loop
        legend_labels = [legend_labels "random seed: " + index * 10]; %#ok<*AGROW>
        mses = [];
        sparsities = [];
        % n_patterns = 1:20;
        for pattern_sp = 50:2:86 %nloc = n_patterns - 1

        disp(rng);
        % CLOSED LOOP PATTERN
        %% PART 1 ------------> Run algorithm with fewer patterns

        % GENERATE HADAMARD CODES, COMPLEMENT INTERLEAVED

        % Define number of orthogonal locations.
        % Usually as n = m - 1 with m a multiple of 4,
        % in any case, m>n patterns are generated, with m as low as possible.
        % The offset defines how orthogonal locations distribute in space; they
        % should be equidistant from each other to minimize scattering crosstalk
        nlocations_and_offset = [(4 * p_index - 1) 1]; % [n offset]
        % nlocations_and_offset(1) = nloc;

        % optional parameter binning to group projection elements
        % binning = [150 150];
        npts = 120 - round(pattern_sp / 2);
        sq_size = 56 - pattern_sp;
        if (sq_size < 1 && pattern_sp < 74)
            sq_size = 36 - round(pattern_sp / 2);
        elseif (sq_size < 1)
            sq_size = 44 - round(pattern_sp / 2);
        end
        
        spacing = sq_size * 2 + round(pattern_sp / 2);

        % image dimensions
        npix_x = 600;
        npix_y = 600;

        % [orig_poisson_patterns, orig_without_complement] = generatePoissonDiscCodesMatrix(nlocations_and_offset, spacing, npts, npix_x, npix_y, sq_size);
        % 
        % figure(1)
        % moviesc(orig_poisson_patterns)
        % title 'illumination patterns (original)'
        % 
        % figure(2)
        % moviesc(orig_without_complement)
        % title 'illumination patterns without interleaved complements (original)'
        % 
        % figure(3)
        % moviesc(vm(sum(orig_without_complement.data, 3)))
        % title 'sum of illumination patterns (original)'

        % evenly spaced poisson disc patterns
        [poisson_patterns, without_complement] = generateEvenlySpacedPoissonDiscCodesMatrix(nlocations_and_offset, spacing, npts, npix_x, npix_y, sq_size);

%         figure(4)
%         moviesc(poisson_patterns)
%         title 'illumination patterns (evenly spaced)'
% 
%         figure(5)
%         moviesc(without_complement)
%         title 'illumination patterns without interleaved complements (evenly spaced)'
% 
%         figure(6)
%         moviesc(vm(sum(without_complement.data, 3)))
%         title 'sum of illumination patterns (evenly spaced)'

        % GENERATE BASE IMAGE (which should have same # pixels as in pattern)
        p = without_complement.data(:, :, 1) > 0;
        numPixels = sum(p, "all");
        [dim1, dim2, dim3] = size(poisson_patterns);

        % centerX1 = 300;
        % centerY1 = 300;
        % radius = round(sqrt(numPixels / pi));
        % centerX2 = 800;
        % centerY2 = 600;
        % radius = 100;

        % centerX = npix_x/2;
        % centerY = npix_y/2;

        centersX = [130, 130, 300, 470, 470];
        centersY = [130, 470, 300, 130, 470];

        distanceImage = generate5BaseImageGradientCircles(centersX, centersY, numPixels, dim1, dim2);
        % distanceImage = generateBaseImageGradientCircle(centerX,centerY,radius, dim1, dim2);

        figure(3);
        d = vm(distanceImage);
        moviesc(d);
        
        %% SAVE FIGURE
        
        fig_name = strcat("Simulations/sparsity/higher_sparsity/origImg_", ...
            int2str(p_index), 'pat_', int2str(index * 10), ...
            "rng_", int2str(pattern_sp), "_spacing");
        saveas(gcf, fig_name+'.png');

        %%

        reconstructedImage = reconstructImage(distanceImage, poisson_patterns, without_complement);

        figure(200);
        moviesc(vm(reconstructedImage));
        title("Image Reconstructed using r patterns");

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
        max_prop = max(props_arr);

        props_arr_temp = props_arr/max_prop;
        % comparison # must be configured depending on image set
        num_clusters = sum(props_arr > .000000000000001);
        obj = bwareafilt(obj, num_clusters);

%         figure(301);
%         moviesc(vm(obj));
%         title('Max Cluster');

        normalized_reconstruction = (reconstructedImage - min(reconstructedImage(:))) / (max(reconstructedImage(:) - min(reconstructedImage(:))));

        [rows, columns, roi_vec] = find(obj .* normalized_reconstruction);


        %% PART 3 ------------> Run algorithm with more patterns using ROI
        % % roi = distanceImage(r_range(1):r_range(2), c_range(1):c_range(2));
        roi = zeros(dim1, dim2);
        for i = 1:length(rows)
            roi(rows(i), columns(i)) = roi_vec(i);
        end
%         figure(202);
%         moviesc(vm(roi));

        % generate hadamard codes
        nlocations_and_offset = [64 14]; % [n offset] % [(4 * index - 1) 1]
        % binning = [43 1]; % optional parameter binning to group projection elements

        [hadamard_patterns, without_complement] = generateEvenlySpacedPoissonDiscCodesVector(nlocations_and_offset, spacing, npts, rows, columns, sq_size);

        reconstructedImage2 = reconstructImageVectorized(roi_vec, hadamard_patterns, without_complement);


        %% PART 4 ------------> Reconstruct to original resolution of image

        [rdim1, rdim2, rdim3] = size(reconstructedImage);
        final = zeros(rdim1, rdim2, rdim3);
        for i = 1:length(rows)
            final(rows(i), columns(i)) = reconstructedImage2(i);
        end


%         figure(204);
%         moviesc(vm(im2double(mat2gray(final))));
%         title("FINAL IMAGE");

        final_with_tvd = reshape(im2double(mat2gray(applyTotalVariationDenoising(final))), dim1, dim2);
        figure(205);
        moviesc(vm(final_with_tvd));
        title("FINAL IMAGE WITH TVD");

        %% MSE, Sparsity Calculation
        
        calc_mse = immse(final_with_tvd, distanceImage);
        disp(calc_mse)
        sparsity = ((npix_x * npix_y) - numPixels) / (npix_x * npix_y);
        disp(sparsity)

        mses(end + 1) = calc_mse;
        sparsities(end + 1) = sparsity;
        
        %% SAVE FIGURE
        
        fig_name = strcat("Simulations/sparsity/higher_sparsity/", ...
            int2str(p_index), 'pat_', int2str(index * 10), ...
            "rng_", num2str(sparsity), "spars");
        saveas(gcf, fig_name+'.png');

        end

        figure(10 + p_index);
        hold on;
        scatter(sparsities, mses);

    end
    figure(10 + p_index);
    hold on;
    legend(legend_labels);
    title("MSE vs Sparsity, Gradient Circles, with Noise, " + (4 * p_index) + " patterns");
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
        image_patterns(:, :, i) = new_img_with_gauss + new_img_with_poiss;
%         image_patterns(:, :, i) = new_img;
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
        image_patterns(:, i) = new_img_with_gauss + new_img_with_poiss;
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

function distanceImage = generateBaseImageCircle(centerX, centerY, radius, dim1, dim2)
    [columnsInImage, rowsInImage] = meshgrid(1:dim2, 1:dim1);

    distanceImage = (rowsInImage - centerY).^2 + (columnsInImage - centerX).^2 <= radius.^2;
    % Normalization
    distanceImage = distanceImage / max(distanceImage(:));
    distanceImage = double(distanceImage) / 10;
end

function distanceImage = generateBaseImageGradientCircle(centerX,centerY,radius, dim1, dim2)
    [columnsInImage, rowsInImage] = meshgrid(1:dim2, 1:dim1);

    circlePixels = (rowsInImage - centerY).^2 + (columnsInImage - centerX).^2 <= radius.^2;
    % Euclidean Distance Transform
    distanceImage = bwdist(~circlePixels);
    % Normalization
    distanceImage = distanceImage / max(distanceImage(:));
    distanceImage = double(distanceImage) / 10;
end

function distanceImage = generate5BaseImageCircles(centerX, centerY, totalPixels, dim1, dim2)
    [columnsInImage, rowsInImage] = meshgrid(1:dim2, 1:dim1);
    
    radius = round(sqrt(totalPixels / (5 * pi)));

    circle1Pixels = (rowsInImage - centerY(1)).^2 + (columnsInImage - centerX(1)).^2 <= radius.^2;
    circle2Pixels = (rowsInImage - centerY(2)).^2 + (columnsInImage - centerX(2)).^2 <= radius.^2;
    circle3Pixels = (rowsInImage - centerY(3)).^2 + (columnsInImage - centerX(3)).^2 <= radius.^2;
    circle4Pixels = (rowsInImage - centerY(4)).^2 + (columnsInImage - centerX(4)).^2 <= radius.^2;
    circle5Pixels = (rowsInImage - centerY(5)).^2 + (columnsInImage - centerX(5)).^2 <= radius.^2;
    distanceImage = circle1Pixels | circle2Pixels | circle3Pixels | circle4Pixels | circle5Pixels;

    % Normalization
    distanceImage = distanceImage / max(distanceImage(:));
    distanceImage = double(distanceImage)/10;
end

function distanceImage = generate5BaseImageGradientCircles(centerX, centerY, totalPixels, dim1, dim2)
    [columnsInImage, rowsInImage] = meshgrid(1:dim2, 1:dim1);
    
    radius = round(sqrt(totalPixels / (5 * pi)));

    circle1Pixels = (rowsInImage - centerY(1)).^2 + (columnsInImage - centerX(1)).^2 <= radius.^2;
    circle2Pixels = (rowsInImage - centerY(2)).^2 + (columnsInImage - centerX(2)).^2 <= radius.^2;
    circle3Pixels = (rowsInImage - centerY(3)).^2 + (columnsInImage - centerX(3)).^2 <= radius.^2;
    circle4Pixels = (rowsInImage - centerY(4)).^2 + (columnsInImage - centerX(4)).^2 <= radius.^2;
    circle5Pixels = (rowsInImage - centerY(5)).^2 + (columnsInImage - centerX(5)).^2 <= radius.^2;
    circlePixels = circle1Pixels | circle2Pixels | circle3Pixels | circle4Pixels | circle5Pixels;
    % Euclidean Distance Transform
    distanceImage = bwdist(~circlePixels);
    % Normalization
    distanceImage = distanceImage / max(distanceImage(:));
    distanceImage = double(distanceImage)/10;
end

function [disc_patterns, without_complement] = generateEvenlySpacedPoissonDiscCodesMatrix(nlocations_and_offset, pt_spacing, npts, npix_x, npix_y, sq_size)
num_frames = nlocations_and_offset(1) + 1;
patterns_logical = false(npix_x, npix_y, num_frames);

numberPts = round((npix_x * npix_y) / ((pt_spacing / 3) ^ 2)) * num_frames;
[all_pts] = round(poissonDisc([npix_x, npix_y], pt_spacing / num_frames, numberPts, 0));

for j=1:num_frames
    indices = randsample(1:length(all_pts), round(length(all_pts) / num_frames), false); 
    pts = all_pts(indices, :);
    
    all_pts(indices, :) = [];

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
end

function [disc_patterns, without_complement] = generateEvenlySpacedPoissonDiscCodesVector(nlocations_and_offset, pt_spacing, npts, index_rows, index_cols, sq_size)
    num_frames = nlocations_and_offset(1) + 1;
    
    roi_rows = ceil((max(index_rows) - min(index_rows)) * 1.15);
    roi_cols = ceil((max(index_cols) - min(index_cols)) * 1.15);
    
    patterns_logical = false(roi_rows, roi_cols, num_frames);
    
    numberPts = round(length(index_rows) / ((pt_spacing / 3) ^ 2)) * num_frames;
    [all_pts] = round(poissonDisc([roi_rows, roi_cols], pt_spacing / num_frames, numberPts, 0));

    for j=1:num_frames+1
        indices = randsample(1:length(all_pts), round(length(all_pts) / num_frames), false); 
        pts = all_pts(indices, :);

        all_pts(indices, :) = [];

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
%     figure(3)
%     moviesc(vm(sum(without_comp.data, 3)))
%     title 'sum of illumination patterns (evenly spaced)'
    
    
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

function [disc_patterns, without_complement] = generatePoissonDiscCodesMatrix(nlocations_and_offset, pt_spacing, npts, npix_x, npix_y, sq_size)
num_frames = nlocations_and_offset(1) + 1;
patterns_logical = false(npix_x, npix_y, num_frames);

for j=1:num_frames
    [pts] = round(poissonDisc([npix_x, npix_y], pt_spacing, npts, 0));

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
%     figure(1)
%     moviesc(patterns)
%     title 'illumination patterns, second round'
% 
%     figure(2)
%     moviesc(without_comp)
%     title 'illumination patterns without interleaved complements, second round'
    
    
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
    img = SB_ATV(reconstructedImg, 0.5);
end