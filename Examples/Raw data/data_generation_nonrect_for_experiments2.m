% Reconstructs a final image using experimental images taken with
% Hadamard optical sectioning microscopy patterns in complement interleaved
% fashion.
%

selftic = tic;

% Dimensions of images taken during experiment
dim1 = 1200;
dim2 = 1200;

%%
% LOAD INTERPOLATED HADAMARD PATTERNS

filename=fullfile('SI/homogenous_randomPattern_r20Pattern_power1_1/Pos0', '*.tif');
directory = dir(filename);

% last image in the directory is not applicable
interpolated_patterns = zeros(dim1, dim2, length(directory)-1);
bin_interpolated_patterns = false(dim1, dim2, length(directory)-1);
for k=1:length(directory)-1
    file = strcat('SI/homogenous_randomPattern_r20Pattern_power1_1/Pos0/', directory(k).name);
    
    I = imread(file);
%     se = strel('disk', 15);
%     filtered = imtophat(I, se);
%     contrastAdjusted = imadjust(filtered);
%     bw = imbinarize(contrastAdjusted);
%     
%     bin_interpolated_patterns(:, :, k) = imbinarize(I(701:1000, 601:900));
      bin_interpolated_patterns(:, :, k) = imbinarize(I, 'adaptive');
%     pattern = mat2gray(imread(file));
%     interpolated_patterns(:, :, k) = I(701:1000, 601:900);
      interpolated_patterns(:, :, k) = I;

end

newInterpolated_pattern=zeros(size(interpolated_patterns));
for i=1:size(interpolated_patterns, 3)
    s = medfilt1(sum(interpolated_patterns(:, :, i), 2),20);
    s_new = repmat(s, [1,size(interpolated_patterns, 2)]);
    temp_pattern=interpolated_patterns(:, :, i) ./ s_new;
    newInterpolated_pattern(:, :, i) = im2double(mat2gray(temp_pattern));
    disp(i)
end

interpolated_patterns = newInterpolated_pattern;

% interp_without_complement = interpolated_patterns(:,:,1:2:end);
% bin_interp_without_complement = bin_interpolated_patterns(:,:,1:2:end);

interp_without_complement = interpolated_patterns(:, :, 1:64);
bin_interp_without_complement = bin_interpolated_patterns(:,:,1:64);

% toint = {interp_without_complement, interpolated_patterns(:, :, 65:128)};
% interpolated_patterns_interleaved = reshape(permute(cell2mat(permute(toint,[1 4 3 2])),[1 2 4 3]),size(toint{1},1),size(toint{1},2),[]);
% interpolated_patterns = interpolated_patterns_interleaved;

% toint = {bin_interp_without_complement, bin_interpolated_patterns(:, :, 65:128)};
% bin_interpolated_patterns_interleaved = reshape(permute(cell2mat(permute(toint,[1 4 3 2])),[1 2 4 3]),size(toint{1},1),size(toint{1},2),[]);
% bin_interpolated_patterns = bin_interpolated_patterns_interleaved;

% figure(2);
% moviesc(vm(interpolated_patterns));
% title('interpolated hadamard patterns, interleaved with complements');

figure(3);
moviesc(vm(bin_interp_without_complement));
title('bin interpolated hadamard patterns without complements');

figure(4);
moviesc(vm(interp_without_complement));
title('interpolated hadamard patterns without complements');

%%
img_name = 'randomPattern_woPhantom_4';
name = strcat('SI/', img_name, '/Pos0');

% LOAD IMAGES THAT WERE TAKEN DURING EXPERIMENT
filename = fullfile(name, '*.tif');
directory = dir(filename);

% last image in the directory is not applicable
image_patterns = zeros(dim1, dim2, length(directory)-1);
for k=1:length(directory)-1
    file = strcat(name, '/', directory(k).name);
%     image = mat2gray(imread(file));
%     image_patterns(:, :, k) = im2double(image);
    image_patterns(:, :, k) = imread(file);
end

toint = {image_patterns(:, :, 1:64), image_patterns(:, :, 65:128)};
image_patterns_interleaved = reshape(permute(cell2mat(permute(toint,[1 4 3 2])),[1 2 4 3]),size(toint{1},1),size(toint{1},2),[]);
image_patterns = image_patterns_interleaved;

figure(1);
moviesc(vm(image_patterns));
title('images taken using interpolated hadamard patterns');

% distanceImage = generateBaseImageGradientCircle(150,150,100, dim1, dim2);

% figure(1);
% d = vm(distanceImage);
% moviesc(d);
% title('base image')

% image_patterns = zeros(dim1, dim2, length(directory)-1);
% for k=1:length(directory)-1
% %     image = mat2gray(imread(file));
% %     image_patterns(:, :, k) = im2double(image);
%     image_patterns(:, :, k) = distanceImage .* bin_interpolated_patterns(:, :,  k);
% %     figure;imagesc(image_patterns(:,:,k));
% end

% figure(2);
% img = vm(image_patterns);
% moviesc(img);
% title('images with hadamard patterns applied')

%%

% % LOAD NON-INTERPOLATED HADAMARD PATTERNS
% load('2560x1600_hadamard_patterns.mat');
% 
% had_patterns = hadamard_pattern_data(1:dim1, 1:dim2, :);
% hadamard_pattern_without_comp = had_patterns(:,:,1:2:end);
% 
% 
% figure(5);
% moviesc(vm(hadamard_pattern_data));
% title('generated hadamard patterns');
% 
% figure(6);
% moviesc(vm(hadamard_pattern_without_comp));
% title('generated hadamard patterns without comp');

%% CLOSED LOOP PATTERN
% PART 1 ------------> Run algorithm with fewer patterns

r = 20;
r_image_patterns = image_patterns(:, :, 1:(r+1));
r_interp_without_complement = bin_interp_without_complement(:, :, 1:(r/2 + 1));

reconstructedImage = reconstructImage(r_image_patterns, r_interp_without_complement);
% reconstructedImage = abs(reconstructedImage);

figure(200);
moviesc(vm(reconstructedImage));
title("Image Reconstructed using r patterns");

% reconstructedImage = im2double(mat2gray(reconstructedImage));
% figure(201);
% moviesc(vm(reconstructedImage));
% title("Image Reconstructed using r patterns, filtered");

% time taken to run first round of estimation
time_to_roi = toc(selftic);

%% SAVE FIGURE & CALC TIME
% fig_name = strcat("Simulations-nonrect_roi_exp/", img_name, '_20_patterns');
% saveas(gcf, fig_name+'.png');
% 
% time_to_roi_3_9umRedbeads_region2_wPhantom_lowcon_3_20_patterns = time_to_roi;
% save('Simulations-nonrect_roi_exp/calc_times', 'time_to_roi_3_9umRedbeads_region2_wPhantom_lowcon_3_20_patterns', '-append');

%%
% PART 2  ------------> Determine ROI using reconstructed image, kmeans

[clusterIndexes, clusterCenters] = kmeans(reconstructedImage(:), 2,...
  'distance', 'sqEuclidean', ...
  'Replicates', 2);
clustering = reshape(clusterIndexes, dim1, dim2);
figure(300);
moviesc(vm(clustering));
title('Kmeans Clustering');

disp(clusterCenters)
[~, indexOfMaxValue] = max(clusterCenters);
obj = clustering == indexOfMaxValue;

props = regionprops(obj, 'FilledArea');
props_arr = cell2mat(struct2cell(props)');
max_prop = max(props_arr);
i = 1;
for j = 1:length(props)
    % comparison # must be configured depending on image set
    if props_arr(i)/max_prop < .000000000000001
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

[rows, columns, roi_vec] = find(obj .* reconstructedImage);

%%
% PART 3 ------------> Run algorithm with more patterns using ROI

interp_without_complement_vec = vectorizeHadamardCodes(interp_without_complement, rows, columns);

total_num_patterns = size(image_patterns, 3);
image_patterns_roi_array = zeros(length(roi_vec), total_num_patterns);
for i = 1:total_num_patterns
    [rows, columns, roi_vec] = find(obj .* image_patterns(:, :, i));
    if length(roi_vec) == length(image_patterns_roi_array(:, i))
        image_patterns_roi_array(:, i) = roi_vec;
    else
        roi_vec(length(image_patterns_roi_array(:, i))) = 0;
        image_patterns_roi_array(:, i) = roi_vec;
    end
end

img_patterns_roi_matrix = zeros(dim1, dim2, size(image_patterns_roi_array, 2));

% view image & illumination patterns
patterns = vm(image_patterns_roi_array);
for d = 1:size(image_patterns_roi_array, 2)
    for i = 1:length(rows)
        img_patterns_roi_matrix(rows(i), columns(i), d) = patterns.data(i, d);
    end
end

figure(20)
moviesc(vm(img_patterns_roi_matrix))
title 'image patterns, second round'

without_complement_matrix = zeros(dim1, dim2, size(interp_without_complement_vec, 2));
for d = 1:size(interp_without_complement_vec, 2)
    for i = 1:length(rows)
        without_complement_matrix(rows(i), columns(i), d) = interp_without_complement_vec(i, d);
    end
end

figure(21)
moviesc(vm(without_complement_matrix))
title 'illumination patterns, second round'

reconstructedImage2 = reconstructImageVectorized(image_patterns_roi_array, interp_without_complement_vec);

%%
% PART 4 ------------> Reconstruct to original resolution of image

[rdim1, rdim2, rdim3] = size(reconstructedImage);
final = zeros(rdim1, rdim2, rdim3);
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

% total time
total_time = toc(selftic);
% time to generate image from roi
time_roi_to_final = total_time - time_to_roi;

%% SAVE FIGURE & CALC TIME
% fig_name = strcat("Simulations-nonrect_roi_exp/", img_name, '_final');
% saveas(gcf, fig_name+'.png');
% 
% time_to_roi_3_9umRedbeads_region2_wPhantom_lowcon_3_final = time_roi_to_final;
% save('Simulations-nonrect_roi_exp/calc_times', 'time_to_roi_3_9umRedbeads_region2_wPhantom_lowcon_3_final', '-append');



function reconstructedImage = reconstructImage(image_patterns, without_complement)
    [imgdim1, imgdim2, ~] = size(image_patterns);
    reconstructedImage = zeros(imgdim1, imgdim2);

    new_image_patterns = vm(image_patterns);
    without_complement = vm(without_complement);
    
%     figure;
%     moviesc(new_image_patterns);
    
%     figure;
%     moviesc(without_complement);
    
    new_image_patterns = new_image_patterns(:, :, 2:end);
    without_complement = without_complement(:, :, 2:end);
    
    figure;
    moviesc(new_image_patterns);
    
    figure;
    moviesc(without_complement);

    % RECONSTRUCT IMAGE
    rootdir = fullfile('Examples/Raw data');
    %%
    res_fpath = fullfile(rootdir,'results.mat');

    %%
    ncomps = 5;
    rng(0,'twister') % ensures multiple runs of eigs() will return repeatable results
    
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
%     w50img = mean(vm(uref*sv',new_image_patterns.imsz)); % widefield image
    h50img = mean(vm(uhad*sv',new_image_patterns.imsz)); % hadamard image

    reconstructedImage(:, :) = h50img;
end

function reconstructedImage = reconstructImageVectorized(image_patterns, without_complement)
    [imgdim1, imgdim2, ~] = size(image_patterns);
    reconstructedImage = zeros(imgdim1, 1);

    new_image_patterns = vm(image_patterns);
    without_complement = vm(without_complement);
    
    %%
    % RECONSTRUCT IMAGE
    rootdir = fullfile('Examples/Raw data');
    res_fpath = fullfile(rootdir,'results.mat');

    %%
    ncomps = 5;
    rng(0,'twister') % ensures multiple runs of eigs() will return repeatable results

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

    reconstructedImage(:, :) = h50img;
end

function hadamard_patterns_vec = vectorizeHadamardCodes(hadamard_patterns_matrix, index_rows, index_cols)
    hdim3 = size(hadamard_patterns_matrix, 3);
    vec_size = length(index_rows);
    
    hadamard_patterns_vec = zeros(vec_size, hdim3);
    
    for d = 1:hdim3
        for i = 1:vec_size
            r = index_rows(i);
            c = index_cols(i);
            hadamard_patterns_vec(i, d) = hadamard_patterns_matrix(r, c, d);
        end
    end
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

function img = applyTotalVariationDenoising(reconstructedImg)
    img = SB_ATV(reconstructedImg, 0.5);
end
