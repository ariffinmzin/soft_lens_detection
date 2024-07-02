
close all;
clear all;
clc;

run = 0; % 0 for single file 1 for all files

if (run==0)
    
    filename = '04261d1016';
% filename = '04261d1142';
% filename = '04261d1036';
% filename = '04261d1033';
% filename = '04851d1283';
% filename = '04851d1336';
% filename = '05044d1163';
% filename = '04397d1854';
% filename = '05015d709'; % no lens
% filename = '05015d795'; % no lens
% filename = '05015d797'; % no lens
% filename = '05156d406';
% filename = '05271d239';
% filename = '05290d130';
% filename = '07013d325'; % error
% filename = '06134d94'; % xde dlm nd_data_training_LG4000

I = imread(strcat(filename,'.tiff'));

%% load data (x y coordinates, radius, etc

load('nd_data_training_LG4000.mat');
n=length(nd_data_training_LG4000);
i=1;
soft_ndiris_training = cell(i,7);
p=length(soft_ndiris_training); 
k=1;

for j=1:n
        
        if (strcmp(nd_data_training_LG4000{j,1},filename)==1)
            
            soft_ndiris_training{k,1} = nd_data_training_LG4000{j,1}; % filename
            soft_ndiris_training{k,2} = nd_data_training_LG4000{j,7}; % pupilx
            soft_ndiris_training{k,3} = nd_data_training_LG4000{j,8}; % pupily
            soft_ndiris_training{k,4} = nd_data_training_LG4000{j,9}; % pupilr
            soft_ndiris_training{k,5} = nd_data_training_LG4000{j,10}; % irisx
            soft_ndiris_training{k,6} = nd_data_training_LG4000{j,11}; % irisy
            soft_ndiris_training{k,7} = nd_data_training_LG4000{j,12}; % irisr
            k=k+1;
            
        end
        
end

%% 

[template, mask, xl_t, yla_t, yla_b, yl_b] = createiristemplate(strcat(soft_ndiris_training{1,1},'.tiff'), soft_ndiris_training(1,2), soft_ndiris_training(1,3), soft_ndiris_training(1,4), soft_ndiris_training(1,5), soft_ndiris_training(1,6), soft_ndiris_training(1,7));

Inoise = imread(strcat('/diagnostics/',soft_ndiris_training{1,1},'.tiff-noise.jpg'));

%# smoothen the image a little with an anisotroic Gaussian
fimg = imfilter(double(Inoise),fspecial('gaussian',[3 1]));
figure, imshow(fimg,[]), title('Gaussian filter');

%# find the lines as local maxima
msk = ones(6);
msk(:,2:5) = 0;
lines = fimg > imdilate(fimg,msk);

%% roi segmentation

imageSize = size(Inoise);
ci = [soft_ndiris_training{1,6},soft_ndiris_training{1,5},soft_ndiris_training{1,7}];     % center and radius of circle ([c_row, c_col, r])
[xx,yy] = ndgrid((1:imageSize(1))-ci(1),(1:imageSize(2))-ci(2));
mask = logical((xx.^2 + yy.^2)>ci(3)^2);
croppedImage = logical(false(size(Inoise)));
croppedImage(:,:,1) = lines(:,:,1).*mask;
figure, imshow(croppedImage), title('inner');

imageSize = size(Inoise);
ci2 = [soft_ndiris_training{1,6},soft_ndiris_training{1,5},soft_ndiris_training{1,7}+40];     % center and radius of circle ([c_row, c_col, r])
[xx,yy] = ndgrid((1:imageSize(1))-ci2(1),(1:imageSize(2))-ci2(2));
mask2 = logical((xx.^2 + yy.^2)<ci2(3)^2);
% mask2 = uint8(mask<ci2(3)^2);
croppedImage2 = logical(false(size(Inoise)));
croppedImage2(:,:,1) = croppedImage(:,:,1).*mask2;
figure, imshow(croppedImage2), title('outer');


%% circular hough transform 

radii = 90:1:200;
h = circle_hough(croppedImage2, radii, 'same', 'normalise');
[peaks, r, c, k] = circle_houghpeaks(h, radii, 'nhoodxy', 15, 'nhoodr', 21, 'npeaks', 1);
% peaks = circle_houghpeaks(h, radii, 'nhoodxy', 297, 'nhoodr', 303, 'npeaks', 2); % ok utk canny matlab



% x1'
% y2'

% figure, imshow(Inoise);
% hold on;
% for peak = peaks
%     [x, y] = circlepoints(peak(3));
%     x+peak(1)
%     y+peak(2)
%     plot(x+peak(1), y+peak(2), 'g-');
% end
% hold off

% yla_t
% peaks(1) - x
% peaks(2) - y
% peaks(3) - radius

[xout_b,yout_b] = linecirc(0,yla_t,peaks(1),peaks(2),peaks(3)); % top intersection with circle
[xout_t,yout_t] = linecirc(0,yla_b,peaks(1),peaks(2),peaks(3)); % bottom intersection with circle

% xout_t
% yout_t; % top intersection
% xout_b
% yout_b; % bottom intersection

% lens_boundary = [x1 y2];
% lens_boundary_cropped_left = [x1 yout_t:yout_b];
% lens_boundary_cropped_right = [x1 yout_t:yout_b];

[x1, y2] = getpeaks(peaks,Inoise);
% 
k1 = find(y2==round(yout_b(1))) % round the number (bottom intersection)
k2 = find(y2==round(yout_t(1))) % round the number (top intersection)

% yla_t % 04261d1016 - 136 (bottom y)
% yla_b % 04261d1016 - 353 (top y)

% [x1, y2] = getpeaks2(peaks,Inoise,yla_t,yla_b,x1);
[x1, y2] = getpeaks2(peaks,Inoise,yla_t,yla_b,x1,soft_ndiris_training{1,1});

% left_soft_lens = k1(1):k2(1); % index 97 - 316 (k1 = 136, k2=356) start
% from top
% a = y2(k1(1)):y2(k2(1))

% index_diff_left = k2(1)-k1(1); % 316 - 97 (substract) array for left
% index_diff_right = k1(2)-k2(2); % 701 - 482

index_diff_left = k1(1)-k2(1); % 316 - 97 array substract for left
index_diff_right = k2(2)-k1(2); % 701 - 482 array substract for right

% [x1, y2] = getpeaks(peaks,Inoise,k1,k2);

% for c = 1:index_diff_left
% left_lens_boundary_y(c) = y2(c+k1(1)-1);
% end

for c = 1:index_diff_left
left_lens_boundary_y(c) = y2(c+k2(1)-1);
end

% for d = 1:index_diff_right
% right_lens_boundary_y(d) = y2(d+k2(2)-1);
% end

for d = 1:index_diff_right
right_lens_boundary_y(d) = y2(d+k1(2)-1);
end

% for e = 1:index_diff_left
% left_lens_boundary_x(e) = x1(e+k1(1)-1);
% end

for e = 1:index_diff_left
left_lens_boundary_x(e) = x1(e+k2(1)-1);
end

% for f = 1:index_diff_right
% right_lens_boundary_x(f) = x1(f+k2(2)-1);
% end

for f = 1:index_diff_right
right_lens_boundary_x(f) = x1(f+k1(2)-1);
end

InoiseA = Inoise(left_lens_boundary_x,left_lens_boundary_y);
% InoiseA(InoiseA==0)=[]; % remove zeros
% InoiseAleft = InoiseA;
% InoiseA = Inoise(1,1);

left_lens_boundary_xy = vertcat(left_lens_boundary_x,left_lens_boundary_y);
right_lens_boundary_xy = vertcat(right_lens_boundary_x,right_lens_boundary_y);



% figure, imshow(Inoise);
% hold on;
% for h = 1:index_diff_right
%     x = left_lens_boundary_x(h);
%     y = left_lens_boundary_y(h);
%   
%     [x1, y2] = getarray(x,y);
%     
%     plot(x1, y2, 'g-');
% end
% hold off
    
% if (y2==round(yout_b(1))) && (x1==round(xout_b(1))) % left boundary
%     
%     
%     
% end
% xout_t =
% 
%   405.9759  240.0241
% 
% 
% yout_t =
% 
%    136   136
% 
% 
% xout_b =
% 
%   419.2912  226.7088
% 
% 
% yout_b =
% 
%    353   353






% out = lineSegmentIntersect([x1, y2],[xl_t, yl_t]);



%% 
    
    
end % run==0


if (run==1)
%% load data (x y coordinates, radius, etc

load('nd_data_training_LG4000.mat');
n=length(nd_data_training_LG4000);
i=1;
soft_ndiris_training = cell(i,7);
p=length(soft_ndiris_training); 
k=1;

for j=1:n
        

            
            soft_ndiris_training{k,1} = nd_data_training_LG4000{j,1}; % filename
            soft_ndiris_training{k,2} = nd_data_training_LG4000{j,7}; % pupilx
            soft_ndiris_training{k,3} = nd_data_training_LG4000{j,8}; % pupily
            soft_ndiris_training{k,4} = nd_data_training_LG4000{j,9}; % pupilr
            soft_ndiris_training{k,5} = nd_data_training_LG4000{j,10}; % irisx
            soft_ndiris_training{k,6} = nd_data_training_LG4000{j,11}; % irisy
            soft_ndiris_training{k,7} = nd_data_training_LG4000{j,12}; % irisr
            k=k+1;
            
%             createiristemplate(strcat(soft_ndiris_training{1,1},'.tiff'), soft_ndiris_training(1,2), soft_ndiris_training(1,3), soft_ndiris_training(1,4), soft_ndiris_training(1,5), soft_ndiris_training(1,6), soft_ndiris_training(1,7));

        
end
%% 

createiristemplate(strcat(soft_ndiris_training{1,1},'.tiff'), soft_ndiris_training(1,2), soft_ndiris_training(1,3), soft_ndiris_training(1,4), soft_ndiris_training(1,5), soft_ndiris_training(1,6), soft_ndiris_training(1,7));

Inoise = imread(strcat(soft_ndiris_training{1,1},'.tiff-noise.jpg'));

%# smoothen the image a little with an anisotroic Gaussian
fimg = imfilter(double(Inoise),fspecial('gaussian',[3 1]));
% figure, imshow(fimg,[]), title('Gaussian filter');

%# find the lines as local maxima
msk = ones(6);
msk(:,2:5) = 0;
lines = fimg > imdilate(fimg,msk);

%% roi segmentation

imageSize = size(Inoise);
ci = [soft_ndiris_training{1,6},soft_ndiris_training{1,5},soft_ndiris_training{1,7}];     % center and radius of circle ([c_row, c_col, r])
[xx,yy] = ndgrid((1:imageSize(1))-ci(1),(1:imageSize(2))-ci(2));
mask = logical((xx.^2 + yy.^2)>ci(3)^2);
croppedImage = logical(false(size(Inoise)));
croppedImage(:,:,1) = lines(:,:,1).*mask;
figure, imshow(croppedImage), title('inner');

imageSize = size(Inoise);
ci2 = [soft_ndiris_training{1,6},soft_ndiris_training{1,5},soft_ndiris_training{1,7}+40];     % center and radius of circle ([c_row, c_col, r])
[xx,yy] = ndgrid((1:imageSize(1))-ci2(1),(1:imageSize(2))-ci2(2));
mask2 = logical((xx.^2 + yy.^2)<ci2(3)^2);
% mask2 = uint8(mask<ci2(3)^2);
croppedImage2 = logical(false(size(Inoise)));
croppedImage2(:,:,1) = croppedImage(:,:,1).*mask2;
figure, imshow(croppedImage2), title('outer');

%% circular hough transform 

radii = 90:1:200;
h = circle_hough(croppedImage2, radii, 'same', 'normalise');
peaks = circle_houghpeaks(h, radii, 'nhoodxy', 15, 'nhoodr', 21, 'npeaks', 1);
% peaks = circle_houghpeaks(h, radii, 'nhoodxy', 297, 'nhoodr', 303, 'npeaks', 2); % ok utk canny matlab

figure, imshow(Inoise);
hold on;
for peak = peaks
    [x, y] = circlepoints(peak(3));
    plot(x+peak(1), y+peak(2), 'g-');
    xy = [x y];
end
hold off


end