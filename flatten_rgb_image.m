function imgGS = flatten_rgb_image(img,plot_flag)
% Convert an MxNxD (color) image array into an MxNx1 grayscale image.
% For example, if image is RGB D=3.
%
% INPUTS
%   img         Original MxNxD image data
%   plot_flag   1=plot image
% OUTPUTS
%   imgGS       Flattened MxNx1 image data (grayscale, GS)

% if number of arguments is less than two, set plot_flag to 0
if nargin<2
    plot_flag = 0; % default to No Plot
end

% if image data is RGB (MxNx3), "flatten" to grayscale (MxNx1)
[M, N, D] = size(img);   % rows , columns , image_depth

if D > 1  % flatten image
    % initialize MxNx1 matrix
    K = zeros(M,N,1);
    % sum all layers (RGB colors) together, then calc. the mean intensity
    for j = 1:D
        % Convert Integer image data to a Double for arithmetic
        K = K + ( double(img(:,:,j)) + 1 );
    end
    
    % calculate the mean intensity at each pixel (divide by D pixels)
    K = K./D;
    
    % convert back into unsigned 8-bit integer format (uint8)
    imgGS = uint8( round(K - 1)); % added ROUND per Matlab suggestion, but think this is redundant
    
else
    % image data was already 1D, just assig to output variable 'imgGS'
    imgGS = img;
end

% show original and flattened (grayscale) image
if plot_flag
    figure;
    subplot(1,2,1); imshow(img);
    title(['original: ',num2str(size(img,1)),'x',num2str(size(img,2)),'x',num2str(size(img,3))])
    subplot(1,2,2); imshow(imgGS);
    title(['flattened: ',num2str(size(imgGS,1)),'x',num2str(size(imgGS,2)),'x',num2str(size(imgGS,3))])
end
