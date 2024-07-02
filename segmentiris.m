% segmentiris - peforms automatic segmentation of the iris region
% from an eye image. Also isolates noise areas such as occluding
% eyelids and eyelashes.
%
% Usage: 
% [circleiris, circlepupil, imagewithnoise] = segmentiris(image)
%
% Arguments:
%	eyeimage		- the input eye image
%	
% Output:
%	circleiris	    - centre coordinates and radius
%			          of the detected iris boundary
%	circlepupil	    - centre coordinates and radius
%			          of the detected pupil boundary
%	imagewithnoise	- original eye image, but with
%			          location of noise marked with
%			          NaN values
%
% Author: 
% Libor Masek
% masekl01@csse.uwa.edu.au
% School of Computer Science & Software Engineering
% The University of Western Australia
% November 2003

function [circleiris, circlepupil, imagewithnoise, xl_t, yla_t, yla_b, yl_b] = segmentiris(eyeimage, pupilx, pupily, pupilr, irisx, irisy, irisr)

% define range of pupil & iris radii

%CASIA
lpupilradius = 28;
upupilradius = 75;
lirisradius = 80;
uirisradius = 150;

%    %LIONS
%    lpupilradius = 32;
%    upupilradius = 85;
%    lirisradius = 145;
%    uirisradius = 169;


% define scaling factor to speed up Hough transform
scaling = 0.4;

reflecthres = 240;

% find the iris boundary
% [row, col, r] = findcircle(eyeimage, lirisradius, uirisradius, scaling, 2, 0.20, 0.19, 1.00, 0.00);
% 
% circleiris = [row col r];
% 
% rowd = double(row);
% cold = double(col);
% rd = double(r);
% 
% irl = round(rowd-rd);
% iru = round(rowd+rd);
% icl = round(cold-rd);
% icu = round(cold+rd);

% modified ---------------------------------------

irisx = cell2mat(irisx);
irisy = cell2mat(irisy);
irisr = cell2mat(irisr);

circleiris = [irisy irisx irisr];

rowd = double(irisy);
cold = double(irisx);
rd = double(irisr);

irl = round(rowd-rd);
iru = round(rowd+rd);
icl = round(cold-rd);
icu = round(cold+rd);

% ------------------------------------------------

imgsize = size(eyeimage);

if irl < 1 
    irl = 1;
end

if icl < 1
    icl = 1;
end

if iru > imgsize(1)
    iru = imgsize(1);
end

if icu > imgsize(2)
    icu = imgsize(2);
end

% to find the inner pupil, use just the region within the previously
% detected iris boundary
imagepupil = eyeimage( irl:iru,icl:icu);

%find pupil boundary
% [rowp, colp, r] = findcircle(imagepupil, lpupilradius, upupilradius ,0.6,2,0.25,0.25,1.00,1.00);
% 
% rowp = double(rowp);
% colp = double(colp);
% r = double(r);
% 
% row = double(irl) + rowp;
% col = double(icl) + colp;
% 
% row = round(row);
% col = round(col);

% modified ---------------------------------------

pupilx = cell2mat(pupilx);
pupily = cell2mat(pupily);
pupilr = cell2mat(pupilr);

rowp = double(pupily); % y-axis
colp = double(pupilx); % x-axis
r = double(pupilr);

row = double(irl) + rowp;
col = double(icl) + colp;

row = round(row);
col = round(col);

% ------------------------------------------------

circlepupil = [row col r];

% set up array for recording noise regions
% noise pixels will have NaN values
imagewithnoise = double(eyeimage);

%find top eyelid
topeyelid = imagepupil(1:(rowp-r),:);
% figure, imshow(topeyelid);
lines = findline(topeyelid);

if size(lines,1) > 0
    [xl yl] = linecoords(lines, size(topeyelid));
    yl_t = double(yl) + irl-1;
    xl_t = double(xl) + icl-1;
    
    yla_t = max(yl_t);
    
    y2 = 1:yla_t; % from 1 to yla (the black region regarded as noise)
    
    ind3 = sub2ind(size(eyeimage),yl_t,xl_t);
    imagewithnoise(ind3) = NaN;
    
    imagewithnoise(y2, xl_t) = NaN;
end

%find bottom eyelid
% bottomeyelid = imagepupil((rowp+r):size(imagepupil,1),:);
bottomeyelid = eyeimage((rowp+r):size(eyeimage,1)-60,icl:icu);
% figure, imshow(bottomeyelid);
lines = findline(bottomeyelid);

display('1');

if size(lines,1) > 0
    
    [xl yl] = linecoords(lines, size(bottomeyelid));
%     yl_b = double(yl)+ irl+rowp+r-2;
    yl_b = double(yl)+rowp+r-2; % modified
    xl_b = double(xl) + icl-1;
    
display('2');   

    yla_b = min(yl_b);
    
    y2 = yla_b:size(eyeimage,1);
    
    ind4 = sub2ind(size(eyeimage),yl_b,xl_b);
    imagewithnoise(ind4) = NaN;
    imagewithnoise(y2, xl_b) = NaN;
    
end

%For CASIA, eliminate eyelashes by thresholding
% ref = eyeimage < 100;
ref = eyeimage < 60; %NDIRIS LG4000
% ref = eyeimage < 30; %NDIRIS AD100
coords = find(ref==1);
imagewithnoise(coords) = NaN;
