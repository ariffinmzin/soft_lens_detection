% add Bag of Words
% for purpose of preparing the featureVector left and right before going to
% the training process
% the second file to run is proposed_method_obj1_280617_BoW_training_runall
% all the features are from training dataset (nd_data_training_LG4000)

close all;
clear all;
clc;

global DIAGPATH
DIAGPATH = 'diagnostics';

% run = 0; % 0 for soft
run = 1; %  1 for no

if (run==0)
    
%% load data (x y coordinates, radius, etc

load('nd_data_training_LG4000.mat');
n=length(nd_data_training_LG4000);
i=1;
soft_ndiris_training = cell(i,7);
p=length(soft_ndiris_training); 
k=1;

for j=1:n
        
try
    
if (strcmp(nd_data_training_LG4000{j,6},'Yes')==1)
            
            soft_ndiris_training{k,1} = nd_data_training_LG4000{j,1}; % filename
            soft_ndiris_training{k,2} = nd_data_training_LG4000{j,7}; % pupilx
            soft_ndiris_training{k,3} = nd_data_training_LG4000{j,8}; % pupily
            soft_ndiris_training{k,4} = nd_data_training_LG4000{j,9}; % pupilr
            soft_ndiris_training{k,5} = nd_data_training_LG4000{j,10}; % irisx
            soft_ndiris_training{k,6} = nd_data_training_LG4000{j,11}; % irisy
            soft_ndiris_training{k,7} = nd_data_training_LG4000{j,12}; % irisr
%             k=k+1;
            
%         end
        


%% 

% imds = imageDatastore(strcat(soft_ndiris_training{k,1},'.tiff-noise.jpg'));

% imds = imageDatastore('04261d1016.tiff-noise.jpg');
% bag = bagOfFeatures(imds);
% img = imread('04261d1016.tiff-noise.jpg');
% img = imread('04261d1016.tiff-noise.jpg');
% featureVector = encode(bag,img);



if exist(strcat(soft_ndiris_training{k,1},'.tiff-lens_segmented.jpg'), 'file') == 0 % if segmented jpg files does not exist, run createiristemplate

display(soft_ndiris_training{k,1});
[template, mask, xl_t, yla_t, yla_b, yl_b] = createiristemplate(strcat(soft_ndiris_training{k,1},'.tiff'), soft_ndiris_training(k,2), soft_ndiris_training(k,3), soft_ndiris_training(k,4), soft_ndiris_training(k,5), soft_ndiris_training(k,6), soft_ndiris_training(k,7));

w = cd;
cd(DIAGPATH);

Inoise = imread(strcat(soft_ndiris_training{k,1},'.tiff-noise.jpg'));

cd(w);

% smoothen the image a little with an anisotroic Gaussian
fimg = imfilter(double(Inoise),fspecial('gaussian',[3 1]));
% figure, imshow(fimg,[]), title('Gaussian filter');

%# find the lines as local maxima
msk = ones(6);
msk(:,2:5) = 0;
lines = fimg > imdilate(fimg,msk);

%% roi segmentation

imageSize = size(Inoise);
% ci = [soft_ndiris_training{1,6},soft_ndiris_training{1,5},soft_ndiris_training{1,7}];     % center and radius of circle ([c_row, c_col, r])
ci = [soft_ndiris_training{k,6},soft_ndiris_training{k,5},soft_ndiris_training{k,7}];     % center and radius of circle ([c_row, c_col, r])
[xx,yy] = ndgrid((1:imageSize(1))-ci(1),(1:imageSize(2))-ci(2));
mask = logical((xx.^2 + yy.^2)>ci(3)^2);
croppedImage = logical(false(size(Inoise)));
croppedImage(:,:,1) = lines(:,:,1).*mask;
% figure, imshow(croppedImage), title('inner');

imageSize = size(Inoise);
% ci2 = [soft_ndiris_training{1,6},soft_ndiris_training{1,5},soft_ndiris_training{1,7}+40];     % center and radius of circle ([c_row, c_col, r])
ci2 = [soft_ndiris_training{k,6},soft_ndiris_training{k,5},soft_ndiris_training{k,7}+40];     % center and radius of circle ([c_row, c_col, r])
[xx,yy] = ndgrid((1:imageSize(1))-ci2(1),(1:imageSize(2))-ci2(2));
mask2 = logical((xx.^2 + yy.^2)<ci2(3)^2);
% mask2 = uint8(mask<ci2(3)^2);
croppedImage2 = logical(false(size(Inoise)));
croppedImage2(:,:,1) = croppedImage(:,:,1).*mask2;
% figure, imshow(croppedImage2), title('outer');


%% circular hough transform 

radii = 90:1:200;
h = circle_hough(croppedImage2, radii, 'same', 'normalise');
[peaks, r, c, k1] = circle_houghpeaks(h, radii, 'nhoodxy', 15, 'nhoodr', 21, 'npeaks', 1); % The first row of PEAKS has the
%   x-coordinates, the second row has the y-coordinates, and the third row
%   has the radii.
% peaks = circle_houghpeaks(h, radii, 'nhoodxy', 297, 'nhoodr', 303, 'npeaks', 2); % ok utk canny matlab

%% lens segmentation for visual and also data for training

[xout_b,yout_b] = linecirc(0,yla_t,peaks(1),peaks(2),peaks(3)); % top intersection with circle
[xout_t,yout_t] = linecirc(0,yla_b,peaks(1),peaks(2),peaks(3)); % bottom intersection with circle
% top n bottom mungkin tersalah tag?

[x1, y2] = getpeaks(peaks,Inoise);
% 
k1 = find(y2==round(yout_b(1))); % round the number (bottom intersection)
k2 = find(y2==round(yout_t(1))); % round the number (top intersection)

% [x1, y2] = getpeaks2(peaks,Inoise,yla_t,yla_b,x1,soft_ndiris_training{1,1});
[x1, y2] = getpeaks2(peaks,Inoise,yla_t,yla_b,x1,soft_ndiris_training{k,1});

index_diff_left = k1(1)-k2(1); % 316 - 97 array substract for left
index_diff_right = k2(2)-k1(2); % 701 - 482 array substract for right

for c = 1:index_diff_left
% left_lens_boundary_y(c) = y2(c+k2(1)-1);
left_lens_boundary_y(c) = x1(c+k2(1)-1);
end

for d = 1:index_diff_right
% right_lens_boundary_y(d) = y2(d+k1(2)-1);
right_lens_boundary_y(d) = x1(d+k1(2)-1);
end

for e = 1:index_diff_left
% left_lens_boundary_x(e) = x1(e+k2(1)-1);
left_lens_boundary_x(e) = y2(e+k2(1)-1);
end

for f = 1:index_diff_right
% right_lens_boundary_x(f) = x1(f+k1(2)-1);
right_lens_boundary_x(f) = y2(f+k1(2)-1);
end

% the data that need to be trained

% left contact lens

w = cd;
cd(DIAGPATH)

% eval(sprintf('Inoise_segmented_left_%s = Inoise(left_lens_boundary_x,left_lens_boundary_y);',soft_ndiris_training{k,1}));
eval(sprintf('Inoise_segmented_left_%s = Inoise(left_lens_boundary_y,left_lens_boundary_x);',soft_ndiris_training{k,1}));
eval(sprintf('Inoise_segmented_left_%s = double(Inoise_segmented_left_%s);',soft_ndiris_training{k,1},soft_ndiris_training{k,1}));
eval(sprintf('Inoise_segmented_left_%s(find(Inoise_segmented_left_%s==0)) = NaN;',soft_ndiris_training{k,1},soft_ndiris_training{k,1}));
% eval(sprintf('Inoise_segmented_left_%s(Inoise_segmented_left_%s==0) = [];',soft_ndiris_training{k,1},soft_ndiris_training{k,1}));
% eval(sprintf('Inoise_segmented_left_%s = Inoise_segmented_left_%s(:);',soft_ndiris_training{k,1},soft_ndiris_training{k,1}));
% eval(sprintf('Inoise_segmented_left_%s = Inoise_segmented_left_%s'';',soft_ndiris_training{k,1},soft_ndiris_training{k,1}));

% eval(sprintf('imds = imageDatastore(''Inoise_segmented_left_%s'');',soft_ndiris_training{k,1}));

eval(sprintf('imwrite(Inoise_segmented_left_%s,''%s_left_lens_region.jpg'')',soft_ndiris_training{k,1},soft_ndiris_training{k,1}));

imds = imageDatastore(strcat(soft_ndiris_training{k,1},'_left_lens_region.jpg')); % lagi elok kalau dpt baca direct Inoise segmented left
bag = bagOfFeatures(imds,'GridStep',[1 1]);
img = readimage(imds, 1);

% bag of words;
eval(sprintf('featureVector_left_%s = encode(bag, img);',soft_ndiris_training{k,1}));

save(strcat('featureVector_left_',soft_ndiris_training{k,1},'.mat'),strcat('featureVector_left_',soft_ndiris_training{k,1}));

% img = imread('04261d1016.tiff-noise.jpg');
% img = imread('04261d1016.tiff-noise.jpg');
% featureVector = encode(bag,img);
% KETUA AKAUNTAN KEMENTERIAN PENDIDIKAN MALAYSIA
% right contact lens

%428489
% eval(sprintf('Inoise_segmented_right_%s = Inoise(right_lens_boundary_x,right_lens_boundary_y);',soft_ndiris_training{k,1}));
eval(sprintf('Inoise_segmented_right_%s = Inoise(right_lens_boundary_y,right_lens_boundary_x);',soft_ndiris_training{k,1}));
eval(sprintf('Inoise_segmented_right_%s = double(Inoise_segmented_right_%s);',soft_ndiris_training{k,1},soft_ndiris_training{k,1}));
eval(sprintf('Inoise_segmented_right_%s(find(Inoise_segmented_right_%s==0)) = NaN;',soft_ndiris_training{k,1},soft_ndiris_training{k,1}));

eval(sprintf('imwrite(Inoise_segmented_left_%s,''%s_right_lens_region.jpg'')',soft_ndiris_training{k,1},soft_ndiris_training{k,1}));

imds = imageDatastore(strcat(soft_ndiris_training{k,1},'_right_lens_region.jpg'));
bag = bagOfFeatures(imds,'GridStep',[1 1]);
img = readimage(imds, 1);

% bag of words
eval(sprintf('featureVector_right_%s = encode(bag, img);',soft_ndiris_training{k,1}));

save(strcat('featureVector_right_',soft_ndiris_training{k,1},'.mat'),strcat('featureVector_right_',soft_ndiris_training{k,1}));

cd(w);

% eval(sprintf('Inoise_segmented_right_%s(Inoise_segmented_right_%s==0) = [];',soft_ndiris_training{k,1},soft_ndiris_training{k,1}));
% eval(sprintf('Inoise_segmented_right_%s = Inoise_segmented_right_%s(:);',soft_ndiris_training{k,1},soft_ndiris_training{k,1}));
% eval(sprintf('Inoise_segmented_right_%s = Inoise_segmented_right_%s'';',soft_ndiris_training{k,1},soft_ndiris_training{k,1}));

% eval(sprintf('left_lens_boundary_xy_%s = horzcat(left_lens_boundary_x,left_lens_boundary_y);', soft_ndiris_training{k,1}));
% eval(sprintf('right_lens_boundary_xy_%s = horzcat(right_lens_boundary_x,right_lens_boundary_y);', soft_ndiris_training{k,1}));
% eval(sprintf('left_right_lens_boundary_xy_%s = horzcat(left_lens_boundary_xy_%s,right_lens_boundary_xy_%s);', soft_ndiris_training{k,1},soft_ndiris_training{k,1},soft_ndiris_training{k,1}));




% save(strcat('left_lens_boundary_xy_',soft_ndiris_training{k,1},'.mat'),strcat('Inoise_segmented_left_',soft_ndiris_training{k,1}));
% save(strcat('right_lens_boundary_xy_',soft_ndiris_training{k,1},'.mat'),strcat('Inoise_segmented_right_',soft_ndiris_training{k,1}));



k=k+1;

end % if (stat == 0)

end % if (strcmp(nd_data_training_LG4000{j,6},'Yes')==1)

catch
    % empty
end

end % for j=1:n  
    
end % run==0




%------------------------------------------------------------------------------------------------------------------------------------------------%


if (run==1) % without lens (added 13/07/17)
    
%% load data (x y coordinates, radius, etc

load('nd_data_training_LG4000.mat');
n=length(nd_data_training_LG4000);
i=1;
no_ndiris_training = cell(i,7);
p=length(no_ndiris_training); 
k=1;

    for j=1:n
        
%         try
    
            if (strcmp(nd_data_training_LG4000{j,6},'No')==1)
            
                no_ndiris_training{k,1} = nd_data_training_LG4000{j,1}; % filename
                no_ndiris_training{k,2} = nd_data_training_LG4000{j,7}; % pupilx
                no_ndiris_training{k,3} = nd_data_training_LG4000{j,8}; % pupily
                no_ndiris_training{k,4} = nd_data_training_LG4000{j,9}; % pupilr
                no_ndiris_training{k,5} = nd_data_training_LG4000{j,10}; % irisx
                no_ndiris_training{k,6} = nd_data_training_LG4000{j,11}; % irisy
                no_ndiris_training{k,7} = nd_data_training_LG4000{j,12}; % irisr
                
                if exist(strcat(no_ndiris_training{k,1},'-lens_segmented.jpg'), 'file') == 0 % if segmented jpg files does not exist, run createiristemplate

                    display(no_ndiris_training{k,1});
%                     [template, mask, xl_t, yla_t, yla_b, yl_b] = createiristemplate(strcat('C:\Users\Arif\Documents\MATLAB\NDIRIS\LG4000\images\',no_ndiris_training{k,1},'.tiff'), no_ndiris_training(k,2), no_ndiris_training(k,3), no_ndiris_training(k,4), no_ndiris_training(k,5), no_ndiris_training(k,6), no_ndiris_training(k,7));
                    [template, mask, xl_t, yla_t, yla_b, yl_b] = createiristemplate(strcat(no_ndiris_training{k,1},'.tiff'), no_ndiris_training(k,2), no_ndiris_training(k,3), no_ndiris_training(k,4), no_ndiris_training(k,5), no_ndiris_training(k,6), no_ndiris_training(k,7));

                    w = cd;
                    cd(DIAGPATH);

                    Inoise = imread(strcat(no_ndiris_training{k,1},'.tiff-noise.jpg'));

                    cd(w);
                    
                    %# smoothen the image a little with an anisotroic Gaussian
                    fimg = imfilter(double(Inoise),fspecial('gaussian',[3 1]));
                    % figure, imshow(fimg,[]), title('Gaussian filter');

                    %# find the lines as local maxima
                    msk = ones(6);
                    msk(:,2:5) = 0;
                    lines = fimg > imdilate(fimg,msk);

                    %% roi segmentation

                    imageSize = size(Inoise);
                    % ci = [soft_ndiris_training{1,6},soft_ndiris_training{1,5},soft_ndiris_training{1,7}];     % center and radius of circle ([c_row, c_col, r])
                    ci = [no_ndiris_training{k,6},no_ndiris_training{k,5},no_ndiris_training{k,7}];     % center and radius of circle ([c_row, c_col, r])
                    [xx,yy] = ndgrid((1:imageSize(1))-ci(1),(1:imageSize(2))-ci(2));
                    mask = logical((xx.^2 + yy.^2)>ci(3)^2);
                    croppedImage = logical(false(size(Inoise)));
                    croppedImage(:,:,1) = lines(:,:,1).*mask;
                    % figure, imshow(croppedImage), title('inner');

                    imageSize = size(Inoise);
                    % ci2 = [soft_ndiris_training{1,6},soft_ndiris_training{1,5},soft_ndiris_training{1,7}+40];     % center and radius of circle ([c_row, c_col, r])
                    ci2 = [no_ndiris_training{k,6},no_ndiris_training{k,5},no_ndiris_training{k,7}+40];     % center and radius of circle ([c_row, c_col, r])
                    [xx,yy] = ndgrid((1:imageSize(1))-ci2(1),(1:imageSize(2))-ci2(2));
                    mask2 = logical((xx.^2 + yy.^2)<ci2(3)^2);
                    % mask2 = uint8(mask<ci2(3)^2);
                    croppedImage2 = logical(false(size(Inoise)));
                    croppedImage2(:,:,1) = croppedImage(:,:,1).*mask2;
                    % figure, imshow(croppedImage2), title('outer');
                    
                    

                    %% circular hough transform 

                    radii = 90:1:200;
                    h = circle_hough(croppedImage2, radii, 'same', 'normalise');
                    [peaks, r, c, k1] = circle_houghpeaks(h, radii, 'nhoodxy', 15, 'nhoodr', 21, 'npeaks', 1);
                    % peaks = circle_houghpeaks(h, radii, 'nhoodxy', 297, 'nhoodr', 303, 'npeaks', 2); % ok utk canny matlab

                    %% lens segmentation for visual and also data for training

                    [xout_b,yout_b] = linecirc(0,yla_t,peaks(1),peaks(2),peaks(3)); % top intersection with circle
                    [xout_t,yout_t] = linecirc(0,yla_b,peaks(1),peaks(2),peaks(3)); % bottom intersection with circle

                    [x1, y2] = getpeaks(peaks,Inoise); % x1,y2 is the segmented lens boundary vector
                    % 
                    k1 = find(y2==round(yout_b(1))); % round the number (bottom intersection)
                    k2 = find(y2==round(yout_t(1))); % round the number (top intersection)

                    % [x1, y2] = getpeaks2(peaks,Inoise,yla_t,yla_b,x1,soft_ndiris_training{1,1});
                    [x1, y2] = getpeaks2(peaks,Inoise,yla_t,yla_b,x1,no_ndiris_training{k,1});

%                     index_diff_left = k1(1)-k2(1); % 316 - 97 array substract for left
%                     index_diff_right = k2(2)-k1(2); % 701 - 482 array substract for right
                    try
                    index_diff_left = k1(1)-k2(1); % 316 - 97 array substract for left
                    catch
                        % empty
                    end
                    
                    try
                    index_diff_right = k2(length(k2))-k1(length(k1)); % 701 - 482 array substract for right
                    % masalah 04267d334 k2 ada 4 elemen [74    75   827
                    % 828]
                    catch
                        % empty
                    end
                    for c = 1:index_diff_left
                        left_lens_boundary_y(c) = x1(c+k2(1)-1);
                    end

                    for d = 1:index_diff_right
                        right_lens_boundary_y(d) = x1(d+k1(2)-1);
                    end

                    for e = 1:index_diff_left
                        left_lens_boundary_x(e) = y2(e+k2(1)-1);
                    end

                    for f = 1:index_diff_right
                        right_lens_boundary_x(f) = y2(f+k1(2)-1);
                    end

                    % the data that need to be trained

                    % left contact lens
                    
                    w = cd;
                    cd(DIAGPATH)
                    
                    try
                        
%                     eval(sprintf('Inoise_segmented_left_%s = Inoise(left_lens_boundary_x,left_lens_boundary_y);',no_ndiris_training{k,1})); % Inoise(index,index)
                    eval(sprintf('Inoise_segmented_left_%s = Inoise(left_lens_boundary_y,left_lens_boundary_x);',no_ndiris_training{k,1})); % Inoise(index,index)
                    eval(sprintf('Inoise_segmented_left_%s = double(Inoise_segmented_left_%s);',no_ndiris_training{k,1},no_ndiris_training{k,1}));
                    eval(sprintf('Inoise_segmented_left_%s(find(Inoise_segmented_left_%s==0)) = NaN;',no_ndiris_training{k,1},no_ndiris_training{k,1}));
                    % eval(sprintf('Inoise_segmented_left_%s(Inoise_segmented_left_%s==0) = [];',soft_ndiris_training{k,1},soft_ndiris_training{k,1}));
                    % eval(sprintf('Inoise_segmented_left_%s = Inoise_segmented_left_%s(:);',soft_ndiris_training{k,1},soft_ndiris_training{k,1}));
                    % eval(sprintf('Inoise_segmented_left_%s = Inoise_segmented_left_%s'';',soft_ndiris_training{k,1},soft_ndiris_training{k,1}));

                    % eval(sprintf('imds = imageDatastore(''Inoise_segmented_left_%s'');',soft_ndiris_training{k,1}));

                    eval(sprintf('imwrite(Inoise_segmented_left_%s,''%s_left_lens_region.jpg'')',no_ndiris_training{k,1},no_ndiris_training{k,1}));

                    imds = imageDatastore(strcat(no_ndiris_training{k,1},'_left_lens_region.jpg'));
                    bag = bagOfFeatures(imds,'GridStep',[1 1]);
                    img = readimage(imds, 1);
                    % featureVector = encode(bag, img);
                    eval(sprintf('featureVector_left_%s = encode(bag, img);',no_ndiris_training{k,1}));

                    save(strcat('featureVector_left_',no_ndiris_training{k,1},'.mat'),strcat('featureVector_left_',no_ndiris_training{k,1}));
                    
                    catch
                        % empty
                    end
                    % img = imread('04261d1016.tiff-noise.jpg');
                    % img = imread('04261d1016.tiff-noise.jpg');
                    % featureVector = encode(bag,img);
                    % KETUA AKAUNTAN KEMENTERIAN PENDIDIKAN MALAYSIA
                    % right contact lens

                    try
                    %428489
                    eval(sprintf('Inoise_segmented_right_%s = Inoise(right_lens_boundary_x,right_lens_boundary_y);',no_ndiris_training{k,1}));
                    eval(sprintf('Inoise_segmented_right_%s = double(Inoise_segmented_right_%s);',no_ndiris_training{k,1},no_ndiris_training{k,1}));
                    eval(sprintf('Inoise_segmented_right_%s(find(Inoise_segmented_right_%s==0)) = NaN;',no_ndiris_training{k,1},no_ndiris_training{k,1}));

                    eval(sprintf('imwrite(Inoise_segmented_left_%s,''%s_right_lens_region.jpg'')',no_ndiris_training{k,1},no_ndiris_training{k,1}));

                    imds = imageDatastore(strcat(no_ndiris_training{k,1},'_right_lens_region.jpg'));
                    bag = bagOfFeatures(imds,'GridStep',[1 1]);
                    img = readimage(imds, 1);
                    % featureVector = encode(bag, img);
                    eval(sprintf('featureVector_right_%s = encode(bag, img);',no_ndiris_training{k,1}));

                    save(strcat('featureVector_right_',no_ndiris_training{k,1},'.mat'),strcat('featureVector_right_',no_ndiris_training{k,1}));

                    catch
                        
                    end
                    cd(w);
                    
                    % eval(sprintf('Inoise_segmented_right_%s(Inoise_segmented_right_%s==0) = [];',soft_ndiris_training{k,1},soft_ndiris_training{k,1}));
                    % eval(sprintf('Inoise_segmented_right_%s = Inoise_segmented_right_%s(:);',soft_ndiris_training{k,1},soft_ndiris_training{k,1}));
                    % eval(sprintf('Inoise_segmented_right_%s = Inoise_segmented_right_%s'';',soft_ndiris_training{k,1},soft_ndiris_training{k,1}));

                    % eval(sprintf('left_lens_boundary_xy_%s = horzcat(left_lens_boundary_x,left_lens_boundary_y);', soft_ndiris_training{k,1}));
                    % eval(sprintf('right_lens_boundary_xy_%s = horzcat(right_lens_boundary_x,right_lens_boundary_y);', soft_ndiris_training{k,1}));
                    % eval(sprintf('left_right_lens_boundary_xy_%s = horzcat(left_lens_boundary_xy_%s,right_lens_boundary_xy_%s);', soft_ndiris_training{k,1},soft_ndiris_training{k,1},soft_ndiris_training{k,1}));

                    % save(strcat('left_lens_boundary_xy_',soft_ndiris_training{k,1},'.mat'),strcat('Inoise_segmented_left_',soft_ndiris_training{k,1}));
                    % save(strcat('right_lens_boundary_xy_',soft_ndiris_training{k,1},'.mat'),strcat('Inoise_segmented_right_',soft_ndiris_training{k,1}));

                   
                    k=k+1;

                end % if (stat == 0)

            end % if (strcmp(nd_data_training_LG4000{j,6},'Yes')==1)

%         catch
            % empty
%         end

    end % for j=1:n  
    
end % run==0



% figure, imshow(Inoise);
% hold on;
% 
% for ab=1:length(right_lens_boundary_x)
% plot(right_lens_boundary_y,right_lens_boundary_x, '-g');
% end
% 
% hold off;
% 
% figure, imshow(Inoise);
% hold on;
% 
% for ab=1:length(left_lens_boundary_x)
% plot(left_lens_boundary_y,left_lens_boundary_x, '-g');
% end
% 
% hold off;
% 
% figure, imshow(Inoise);
% hold on;
% 
% for ab=1:length(left_lens_boundary_x)
% plot(left_lens_boundary_x,left_lens_boundary_y, '-g');
% end
% 
% hold off;
