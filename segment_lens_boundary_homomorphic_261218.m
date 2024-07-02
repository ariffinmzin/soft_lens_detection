% created on 261218 to segment the soft lens boundary from enhanced images of HE, CLAHE and
% homomorphic filter resulted from apply_image_enhancements_261218.m

clear all
close all
clc

% list images
load('nd_data_training_LG4000.mat');
n=length(nd_data_training_LG4000);
i=1;
soft_ndiris_training = cell(i,7);
p=length(soft_ndiris_training); 
k=1;

for j=1:n

    if (strcmp(nd_data_training_LG4000{j,6},'Yes')==1)
            
    soft_ndiris_training{k,1} = nd_data_training_LG4000{j,1}; % filename
    soft_ndiris_training{k,2} = nd_data_training_LG4000{j,7}; % pupilx
    soft_ndiris_training{k,3} = nd_data_training_LG4000{j,8}; % pupily
    soft_ndiris_training{k,4} = nd_data_training_LG4000{j,9}; % pupilr
    soft_ndiris_training{k,5} = nd_data_training_LG4000{j,10}; % irisx
    soft_ndiris_training{k,6} = nd_data_training_LG4000{j,11}; % irisy
    soft_ndiris_training{k,7} = nd_data_training_LG4000{j,12}; % irisr
    
    if exist(strcat(soft_ndiris_training{k,1},'.tiff-homomorphic-lens_segmented.jpg'), 'file') == 0 % if segmented jpg files does not exist, run code below
        
        display(soft_ndiris_training{k,1}); 
        
        I = imread(strcat(soft_ndiris_training{k,1},'.tiff-homomorphic.jpg')); % homomorphic image
        I_original = imread(strcat(soft_ndiris_training{k,1},'.tiff')); % original image
        
        %% smoothen the image a little with an anisotroic Gaussian

        fimg = imfilter(double(I),fspecial('gaussian',[3 1]));
        % figure, imshow(fimg,[]), title('Gaussian filter');

        %% find the lines as local maxima

        msk = ones(6);
        msk(:,2:5) = 0;
        lines = fimg > imdilate(fimg,msk);

        %% roi segmentation

        imageSize = size(I);
        % ci = [soft_ndiris_training{1,6},soft_ndiris_training{1,5},soft_ndiris_training{1,7}];     % center and radius of circle ([c_row, c_col, r])
        ci = [soft_ndiris_training{k,6},soft_ndiris_training{k,5},soft_ndiris_training{k,7}];     % center and radius of circle ([c_row, c_col, r])
        [xx,yy] = ndgrid((1:imageSize(1))-ci(1),(1:imageSize(2))-ci(2));
        mask = logical((xx.^2 + yy.^2)>ci(3)^2);
        croppedImage = logical(false(size(I)));
        croppedImage(:,:,1) = lines(:,:,1).*mask;
        % figure, imshow(croppedImage), title('inner');

        imageSize = size(I);
        % ci2 = [soft_ndiris_training{1,6},soft_ndiris_training{1,5},soft_ndiris_training{1,7}+40];     % center and radius of circle ([c_row, c_col, r])
        ci2 = [soft_ndiris_training{k,6},soft_ndiris_training{k,5},soft_ndiris_training{k,7}+40];     % center and radius of circle ([c_row, c_col, r]) add another 40 pixels
        [xx,yy] = ndgrid((1:imageSize(1))-ci2(1),(1:imageSize(2))-ci2(2));
        mask2 = logical((xx.^2 + yy.^2)<ci2(3)^2);
        % mask2 = uint8(mask<ci2(3)^2);
        croppedImage2 = logical(false(size(I)));
        croppedImage2(:,:,1) = croppedImage(:,:,1).*mask2;
        % figure, imshow(croppedImage2), title('outer');
    
        %% circular hough transform
        
        radii = 90:1:200;
        h = circle_hough(croppedImage2, radii, 'same', 'normalise');
        [peaks, r, c, k1] = circle_houghpeaks(h, radii, 'nhoodxy', 15, 'nhoodr', 21, 'npeaks', 1); % The first row of PEAKS has the
        %   x-coordinates, the second row has the y-coordinates, and the third row
        %   has the radii.
        % peaks = circle_houghpeaks(h, radii, 'nhoodxy', 297, 'nhoodr', 303, 'npeaks', 2); % ok utk canny matlab
        
        %% lens boundary segmentation
        
        [xout_b,yout_b] = linecirc(0,yla_t,peaks(1),peaks(2),peaks(3)); % top intersection with circle
        [xout_t,yout_t] = linecirc(0,yla_b,peaks(1),peaks(2),peaks(3)); % bottom intersection with circle
        
        [x1, y2] = getpeaks3(peaks,I_original,soft_ndiris_training{k,1}); % own function
        
        k=k+1;
        
    end
    
    end
        
end