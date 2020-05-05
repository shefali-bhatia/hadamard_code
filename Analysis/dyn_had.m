% Copyright 2018 Vicente Parot
% 
% Permission is hereby granted, free of charge, to any person obtaining a
% copy of this software and associated documentation files (the
% "Software"), to deal in the Software without restriction, including
% without limitation the rights to use, copy, modify, merge, publish,
% distribute, sublicense, and/or sell copies of the Software, and to permit
% persons to whom the Software is furnished to do so, subject to the
% following conditions:      
% 
% The above copyright notice and this permission notice shall be included
% in all copies or substantial portions of the Software.    
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
% MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
% NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
% DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
% OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
% USE OR OTHER DEALINGS IN THE SOFTWARE.      
%
function [sv, uhad, uref, uu] = dyn_had(obj_data, cal_data, ncomps)
%dyn_had   Dynamic hadamard demodulation. 
%   dyn_had(obj_data, cal_data, ncomps) analyzes a movie of a dynamic
%   object illuminated with periodic interleaved patterned illumination
%   using hadamard codes. Then, it calculates a low rank estimate of the
%   dynamic object movie as if illuminated uniformly, and a
%   hadamard-demodulated  optical section of the low rank dynamic object.
%
%   Input:
%       obj_data:   Matrix of size [npix nframes] with the patterned
%                   illumination movie, where npix is the number of pixels
%                   in the movie and nframes the number of timepoints. Odd
%                   frames should have been illuminated by a periodic
%                   repetition of the patterns in cal_data, while even
%                   frames should have been illuminated by their binary
%                   complement.
%                   EACH COL = pixels of image illuminated in frame
%       cal_data:	Matrix of size [npix hlen] with the calibration
%                   dataset, where hlen is the length of the hadamard code.
%       ncomps:     Number of PCA components to represent the dataset. This
%                   is an upper bound on the rank of the movie, so it
%                   should be enough to include as many as distinct
%                   dynamical objects are represented in the movie.
%
%   Output: data follows the notation of pca_eig, where a matrix M is
%       approximated by a truncated singular value decomposition U*S*V'.
%       sv:     Temporal functions of the dataset. The same temporal
%               functions are common to all output movies. These column
%               vectors incorporate their eigenvalue and thus sv is
%               orthogonal but not orthonormal. 
%       uhad:   Spatial functions demodulated from the patterned dataset,
%               representing the optical section. uhad*sv' is a matrix
%               of size [npix nframes] with an estimated movie of the
%               uniformly illuminated dynamical object optical section.
%       uref:   Spatial functions averaged from the patterned dataset,
%               representing the widefield reference. uref*sv' is a matrix
%               of size [npix nframes] with an estimated movie of the
%               uniformly illuminated dynamical object with no sectioning.
%       uu:     Spatial functions calculated by pca_eig from a filtered
%               version of obj_data. If correctly estimated, uref should be
%               similar to uu. uu*sv' approximates a filtered version of
%               obj_data that serves as a uniform illuminated reference for
%               extraction of dynamical signals.
%
%   See also vm, pca_eig, hadamard_bincode.
% 
%   2018 Vicente Parot
%   Cohen Lab - Harvard University

    %%

    mov_obj = vm(obj_data,[size(obj_data,1) 1]); % recast input matrix as vectorized movie with image dimensions [npix 1] (otherwise would need to know image dimensions, or require object inputs)
    hlen = size(cal_data,2); % hadamard code length
    interlen = hlen*2; % interleaved code length, double of hadamard length

    %% Filtered uniform movie
    % filter removes residual periodic illumination artifact from imperfect
    % illumination pattern projection
    ker = [0.5; ones([hlen-1 1]); 0.5]; % filter kernel, sliding window average
    ker = permute(ker,[3 2 1])./sum(ker);
    mov_obj_pairs = mov_obj.blnfun(@mean,2); % widefield movie from sum of pairs
    lf = mov_obj_pairs.imfilter(ker,'replicate'); % movie with low temporal frequency dynamics
    lf(lf<2) = 2; % avoid division errors
    nk = evnfun(mov_obj_pairs,@mean,hlen)./evnfun(lf,@mean,hlen); % extract normalized periodic component from movie
    nk = nk([hlen/2+1:end 1:hlen/2]); % rearrange frames in appropriate order
    % correct high-speed periodic intensity fluctuations. input movie must
    % have length of integer multiple of the pattern sequence. mov_uni has
    % half the frames of the input movie.
%     mov_uni = mov_obj_pairs./nk.repmat([1 1 mov_obj_pairs.frames/hlen]);
    [~, ~, lfdim3] = size(lf);
    mov_uni = mov_obj(:,:,1:lfdim3);
    % mov_uni should not have fluctuations that are periodic with the
    % period of the pattern sequence 
    
    %% Autocorrelation Processing
    
    %%
    % Extract PSF function from frames
    data = mov_uni(:, :, :).data;

    reshaped_frames = reshape(data, 600, 600, []);
    
%     for i = 1:size(reshaped_frames, 3)
    i = 1;
    frame = reshaped_frames(:, :, i);
%     frame_norm = (frame(:) - min(frame(:))) / max((frame(:) - min(frame(:))));
%     frame_norm = reshape(frame_norm, 600, 600);
    
%     counts = imhist(frame, 100);
%     
%     T = otsuthresh(counts);
%     thresh = imbinarize(frame, T*3.5);
%     figure();subplot(1,2,1);imagesc(frame);
%     subplot(1,2,2);imagesc(thresh);
    
%     frame_bin = imbinarize(frame, 0.1);
%     
%     figure; 
%     moviesc(vm(frame_bin));
%     title("binarized frame");
   
    [idx, ~] = kmeans(frame(:), 2, 'distance', 'cityblock', ...
        'Replicates', 2, 'Start', 'cluster');

%     [idx, ~] = dbscan(frame(:), 3, 5, 'Distance', 'squaredeuclidean');
    
    clustering = reshape(idx, 600, 600);
    figure(300);
    moviesc(vm(clustering));
    title("clustering");
    
    B = bwboundaries(clustering, 'noholes');

    numClusters = size(B, 1);
    coord = zeros(numClusters, 2);
    for j = 1:numClusters
        coord(j, 1) = ceil(mean(B{j}(:, 1)));
        coord(j, 2) = ceil(mean(B{j}(:, 2)));
    end
% end

    r = 10;
    
%     [cx, cy] = ind2sub(size(thresh), find(thresh));
%     [cyNum,cyVal] = hist(cy, round(numel(cy)/10));
%     [~, Indtemp] = sort(cyNum,'descend');
%     temp2 = cyVal(Indtemp(1:4));
    
    idx_temp = zeros(r*2 + 1, 2);
    for i = 1:size(coord, 1)            %temp2
        % get radius box around determined points
        idx_temp(:, 1) = [coord(i, 1)-r:coord(i, 1)+r];
        idx_temp(:, 2) = [coord(i, 2)-r:coord(i, 2)+r];
%         idx_temp = find(cy>(temp2(i)-10) & cy<(temp2(i)+10));
        H = zeros(size(idx_temp, 1), r*2 + 1);
        for j = 1:size(idx_temp, 1)         % numel(idx_temp)
            start = idx_temp(j, 1);
            end_ = idx_temp(j, 2)+r;
            if start > 600
                start = 600;
            end
            if end_ > 600
                end_ = 600;
            end
            % check where points are filled in
            a = clustering(start, idx_temp(j, 2)-r:end_);
            [~, temp] = max(a);
            
            start_y = idx_temp(j, 2)-r+(temp-r);
            end_y = idx_temp(j, 2)+temp;
            if start_y < 1
                temp2 = clustering(start, 1:idx_temp(j, 2)+temp);
                to_append = zeros(1 - start_y, 1);
                  
                H(j, :) = cat(2, temp2, to_append');
            end
            if end_y > 600
                temp2 = clustering(start, start_y:600);
                to_append = ones(end_y - 600, 1) * 600;
                  
                H(j, :) = cat(2, temp2, to_append');
            end
            if start_y > 1 && end_y < 600
                H(j, :) = clustering(start, start_y:end_y);
            end

            
%             a = thresh(cx(idx_temp(j)), cy(idx_temp(j))-r:cy(idx_temp(j))+r);
%             [~,temp3]=max(a);
%             H(j, :) = thresh(cx(idx_temp(j)),cy(idx_temp(j))-r+(temp3-r):cy(idx_temp(j))+r+(temp3-r));
%             
            % update psf
            H(j, :) = (H(j,:) - min(H(j,:))) / max(H(j,:) - min(H(j,:)));
        end
        
        if i == 1
            H1 = H;
        else    
            H1 = cat(1, H1, H);
        end
    end
    
    H1(~any(~isnan(H1), 2),:)=[];
    
    xforH = repmat(-r:r, [size(H1, 1), 1]);
    
    [x,y] = prepareCurveData(xforH, H1);
    
    % second order gaussian fit (estimate scattering)
    [fitobject, ~] = fit(x, y, 'gauss2');
    
    xPlot = linspace(x(1), x(end), 100);
    figure;plot(xPlot, fitobject(xPlot));
    
    deconvFrame = deconvlucy(frame, fitobject(xPlot));
    figure(301);moviesc(vm(deconvFrame));title("deconvoluted frame");
    

%%
    corr = xcorr(data);
    psf = fspecial('gaussian', 600, 10);
%     psf = imgaussfilt3(reshaped_frames);
%     estimated_nsr = 0.0001 / var(corr(:));
%     temp = deconvwnr(reshaped_frames, psf, estimated_nsr);


    [temp, est_psf] = deconvblind(reshaped_frames, psf);
    
    figure; moviesc(vm(est_psf)); title("estimated psf function");
    figure; moviesc(vm(temp)); title("deconv img");
    
    
    mov_uni_autocorr = reshape(temp, 360000, 1, []);
    
    % Iterative phase retrieval (Gerchberg-Saxton algorithm)
%     input_intensity = mov_uni(:, :, :).data;
    A = fftshift(ifft2(fftshift(data)));
    [mdim1, mdim2, mdim3] = size(mov_uni_autocorr);
    result = zeros(mdim1, mdim2, mdim3);
    for i=1:25
      B = abs(mov_uni_autocorr) .* exp(1i*angle(A));
      C = fftshift(fft2(fftshift(B)));
      D = abs(mov_uni(:, :, :).data) .* exp(1i*angle(C));
      A = fftshift(ifft2(fftshift(D)));
      result = abs(C);
    end
    
    result_reshaped = reshape(result, 600, 600, []);
    figure; moviesc(vm(result_reshaped)); 
    title("img after phase retrieval alg");
    
    mov_uni = vm(result);
    
    %% Direct full estimation
    idx_map = mod((1:mov_obj.frames)-1,interlen)+1; % maps the frame index of the calibration data into the movie
    selftic = tic; % mark time before pca
    [uu, su, vu] = svds(mov_uni(:,:) - mean(mov_uni(:,:)),ncomps); % mov_uni is a uniform movie previously estimated.
    tpca = toc(selftic); % time after factorization
    comp_uu = uu;
    uu = real(uu); % cast to avoid complex datatype from complex eigenvalues 
    sv = real(vu*su); % mov_uni(:,:) === uu*sv' where uu columns are orthonormal
    uhad = zeros([prod(mov_obj.imsz) ncomps]); % allocate memory for results
    uref = zeros([prod(mov_obj.imsz) ncomps]);
    residual_mov = mov_obj; % initialize residual as input movie
    disp 'estimating components ...'
    
    [coeff, score, latent] = pca(mov_uni(:,:), 'Algorithm', 'svd');
    
%     reshaped_frames = reshape(mov_uni(:, :, :).data, 600, 600, []);
%     corr = xcorr(mov_uni(:, :, :).data);
%     psf = fspecial('disk', 10);
%     estimated_nsr = 0.0001 / var(corr(:));
%     temp = deconvwnr(reshaped_frames, psf, estimated_nsr);
    
    pc = comp_uu * su;
    for comp = 1:ncomps % for each component of the USV decomposition
        sv3d = reshape(sv([1:end; 1:end],comp),1,1,[]); % arrange time trace along dimension 3
        ui = evnfun(residual_mov.*sv3d,@sum,interlen)./evnfun((mov_uni(1,1,[1:end; 1:end])*0+1).*sv3d.^2,@sum,interlen); % interpolate
       
          
        [dim1, dim2, dim3] = size(ui);
       
%         figure(99+comp);
%         imhist(ui_reshaped);
%         title(comp);
        
%         for i = 1:dim3
%             current_frame = data(:, :, i);
%             
%             if mean(current_frame(:)) > 0
%                 disp(mean(current_frame(:)));
%                 current_frame = abs(current_frame);
%                 thresh = graythresh(current_frame);
%                 disp(thresh);
%                 curr_frame_thresh = wthresh(current_frame, 's', thresh);
%                 data(:, :, i) = curr_frame_thresh;
%             elseif mean(current_frame(:)) < 0
%                 disp(mean(current_frame(:)));
%                 current_frame = abs(current_frame);
%                 thresh = graythresh(current_frame);
%                 disp(thresh);
%                 curr_frame_thresh = wthresh(current_frame, 's', thresh);
%                 data(:, :, i) = -1 .* curr_frame_thresh;
%             end
%         end
        
        ui_reshaped = reshape(ui.data, 600, 600, dim3);
        figure(9+comp);
        moviesc(vm(ui_reshaped));
        title(comp);
        
        
        uhad(:,comp) = mean((ui(:,1:2:end) - ui(:,2:2:end)).*cal_data,2); % hadamard demodulation
        uref(:,comp) = mean(ui(:,:),2); % widefield reference
        residual_mov = residual_mov - ui(idx_map).*sv3d; % update residual
    end
    t_all = toc(selftic); % total time
    t_est = t_all-tpca; % time for estimation only 
    disp(['estimation took ' num2str(t_est) ' s, total ' num2str(t_all) ' s']);
