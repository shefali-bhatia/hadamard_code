% Reconstructs a final image using experimental images taken with
% Hadamard optical sectioning microscopy patterns in complement interleaved
% fashion.
%

selftic = tic;

% Dimensions of images taken during experiment
dim1 = 600;
dim2 = 600;

%%
% LOAD INTERPOLATED HADAMARD PATTERNS
filename=fullfile('SI/homogeneousSample_1/Pos0/', '*.tif');
directory = dir(filename);

% last image in the directory is not applicable
interpolated_patterns = zeros(dim1, dim2, length(directory)-1);
bin_interpolated_patterns = false(dim1, dim2, length(directory)-1);
for k=1:length(directory)-1
    file = strcat('SI/homogeneousSample_1/Pos0/', directory(k).name);
    I = imread(file);
    bin_interpolated_patterns(:, :, k) = imbinarize(I, 'adaptive');
    interpolated_patterns(:, :, k) = I;
end

newInterpolated_pattern=zeros(size(interpolated_patterns));
for i=1:size(interpolated_patterns, 3)
    s = medfilt1(sum(interpolated_patterns(:, :, i), 2),20);
    s_new = repmat(s, [1,size(interpolated_patterns, 2)]);
    temp_pattern=interpolated_patterns(:, :, i) ./ s_new;
    newInterpolated_pattern(:, :, i) = im2double(mat2gray(temp_pattern));
%     disp(i)
end

interpolated_patterns = newInterpolated_pattern;
interp_without_complement = interpolated_patterns(:, :, 1:64);
bin_interp_without_complement = bin_interpolated_patterns(:,:,1:64);

% figure(3);
% moviesc(vm(bin_interp_without_complement));
% title('bin interpolated hadamard patterns without complements');

% figure(4);
% moviesc(vm(interp_without_complement));
% title('interpolated hadamard patterns without complements');

%%

original = double(imread('SI/3.9umRedbeads_region2_1/Pos0/img_000000000_Default_000.tif'));
figure;moviesc(vm(original));

img_name = '3.9umRedbeads_region2_1';
name = strcat('SI/', img_name, '/Pos0');

% LOAD IMAGES THAT WERE TAKEN DURING EXPERIMENT
filename = fullfile(name, '*.tif');
directory = dir(filename);

% last image in the directory is not applicable
image_patterns = zeros(dim1, dim2, length(directory)-1);
snr = zeros(length(directory)-1, 1);
for k=1:length(directory)-1
    file = strcat(name, '/', directory(k).name);
    
    orig_img = imread(file);
    
    % Add noise (for simulations)
    new_img_with_gauss = imnoise(orig_img,'gaussian');
%     new_img_with_poiss = imnoise(orig_img, 'poisson');
%     new_img = new_img_with_gauss + new_img_with_poiss;
    
    % Normalize and filter out background of images
    norm_img = im2double(mat2gray(orig_img));
    thresh = graythresh(norm_img);
    filtered_img = wthresh(norm_img, 's', thresh);

    image_patterns(:, :, k) =  orig_img;
end

toint = {image_patterns(:, :, 1:64), image_patterns(:, :, 65:128)};
image_patterns_interleaved = reshape(permute(cell2mat(permute(toint,[1 4 3 2])),[1 2 4 3]),size(toint{1},1),size(toint{1},2),[]);
image_patterns = image_patterns_interleaved;

% figure(1);
% moviesc(vm(image_patterns));
% title('images taken using interpolated hadamard patterns');

%% CLOSED LOOP PATTERN
% PART 1 ------------> Run algorithm with fewer patterns

r = 20;
r_image_patterns = image_patterns(:, :, 1:(r+1));
r_interp_without_complement = bin_interp_without_complement(:, :, 1:(r/2 + 1));

reconstructedImage = reconstructImage(r_image_patterns, r_interp_without_complement);

% figure(200);
% moviesc(vm(reconstructedImage));
% title("Image Reconstructed using r patterns");

% time taken to run first round of estimation
time_to_roi = toc(selftic);

%% SAVE FIGURE & CALC TIME
% fig_name = strcat("Simulations-nonrect_roi_exp/", img_name, '_20_patterns');
% saveas(gcf, fig_name+'.png');
% 
% time_to_roi_3_9umRedbeads_region2_wPhantom_lowcon_3_20_patterns = time_to_roi;
% save('Simulations-nonrect_roi_exp/calc_times', 'time_to_roi_3_9umRedbeads_region2_wPhantom_lowcon_3_20_patterns', '-append');

time_20 = time_to_roi;
mse_20 = immse(original, reconstructedImage);

%%
% PART 2  ------------> Determine ROI using reconstructed image, kmeans

[clusterIndexes, clusterCenters] = kmeans(reconstructedImage(:), 2,...
  'distance', 'sqEuclidean', ...
  'Replicates', 2);
clustering = reshape(clusterIndexes, dim1, dim2);
% figure(300);
% moviesc(vm(clustering));
% title('Kmeans Clustering');

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

% figure(301);
% moviesc(vm(obj));
% title('Max Cluster');

normalized_reconstruction = (reconstructedImage - min(reconstructedImage(:))) / (max(reconstructedImage(:) - min(reconstructedImage(:))));

[rows, columns, roi_vec] = find(obj .* normalized_reconstruction);

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

reconstructedImage2 = reconstructImageVectorized(image_patterns_roi_array, interp_without_complement_vec);

%%
% PART 4 ------------> Reconstruct to original resolution of image

[rdim1, rdim2, rdim3] = size(reconstructedImage);
final = zeros(rdim1, rdim2, rdim3);
for i = 1:length(rows)
    final(rows(i), columns(i)) = reconstructedImage2(i);
end


final_norm = im2double(mat2gray(final));


% figure(204);
% moviesc(vm(final_norm));
% title("FINAL IMAGE");

final_with_tvd = reshape(applyTotalVariationDenoising(final), dim1, dim2);
final_with_tvd_norm = im2double(mat2gray(final_with_tvd));
% figure(205);
% moviesc(vm(final_with_tvd_norm));
% title("FINAL IMAGE WITH TVD AND DECONVOLUTION");


% % figure(206);
% % filtered_background = im2double(mat2gray(im2double(mat2gray(image_patterns)) + final_with_tvd));
% % moviesc(vm(filtered_background > 0.5));
% % title("Final Reconstruction with TVD and Deconvolution");
% % 
% % 
% % figure(206);
% % filtered_background2 = im2double(mat2gray(im2double(mat2gray(image_patterns)) + final));
% % moviesc(vm(filtered_background > 0.5));
% % title("Final Reconstruction with only SVD");

% temp = im2double(mat2gray(image_patterns)) + final_with_tvd_norm;
% figure;moviesc(vm(wthresh(temp, 's', graythresh(temp))));
% 
% temp = im2double(mat2gray(image_patterns)) + final_norm;
% figure;moviesc(vm(wthresh(temp, 's', graythresh(temp))));

% figure;moviesc(vm(wthresh(final_with_tvd_norm, 's', graythresh(final_with_tvd_norm))));title("FINAL IMAGE WITH TVD AND DECONVOLUTION");

reg = wthresh(final_norm, 's', graythresh(final_norm));
wtvd = wthresh(final_with_tvd_norm, 's', graythresh(final_with_tvd_norm));
for i=1:3
    wtvd = wthresh(wtvd, 's', graythresh(wtvd));
end

figure;moviesc(vm(im2double(mat2gray(wthresh(reg, 's', graythresh(reg))))));title("Final Reconstruction with only SVD");

figure;moviesc(vm(im2double(mat2gray(wtvd))));title("Final Reconstruction with TVD and Deconvolution");

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
    
    new_image_patterns = new_image_patterns(:, :, 2:end);
    without_complement = without_complement(:, :, 2:end);
    
%     figure;
%     moviesc(new_image_patterns);
%     
%     figure;
%     moviesc(without_complement);

    % RECONSTRUCT IMAGE
    rootdir = fullfile('Examples/Raw data');
    %%
    res_fpath = fullfile(rootdir,'results.mat');

    %%
    ncomps = 1;
    rng(0,'twister') % ensures multiple runs of eigs() will return repeatable results
    
    [sv, uhad, uref, uu] = dyn_had(new_image_patterns(:, :), without_complement(:,:), ncomps);
    %%
    save(res_fpath,...
        'sv', ...
        'uhad', ...
        'uref', ...
        'uu', ...
        'ncomps')
    %%
    % extract time-integrated images 
    h50img = mean(vm(uhad*sv', new_image_patterns.imsz)); % hadamard image

    reconstructedImage(:, :) = h50img;
end

function reconstructedImage = reconstructImageVectorized(image_patterns, without_complement)
    [imgdim1, ~, ~] = size(image_patterns);
    reconstructedImage = zeros(imgdim1, 1);

    new_image_patterns = vm(image_patterns);
    without_complement = vm(without_complement);
    
    %%
    % RECONSTRUCT IMAGE
    rootdir = fullfile('Examples/Raw data');
    res_fpath = fullfile(rootdir,'results.mat');

    %%
    ncomps = 1;
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
    % extract time-integrated images 
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

function img = applyTotalVariationDenoising(reconstructedImg)
    img = SB_ATV(reconstructedImg, 0.5);
end
