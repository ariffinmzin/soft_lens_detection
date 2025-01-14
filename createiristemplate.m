% createiristemplate - generates a biometric template from an iris in
% an eye image.
%
% Usage: 
% [template, mask] = createiristemplate(eyeimage_filename)
%
% Arguments:
%	eyeimage_filename   - the file name of the eye image
%
% Output:
%	template		    - the binary iris biometric template
%	mask			    - the binary iris noise mask
%
% Author: 
% Libor Masek
% masekl01@csse.uwa.edu.au
% School of Computer Science & Software Engineering
% The University of Western Australia
% November 2003

% function [template, mask] = createiristemplate(eyeimage_filename)
function [template, mask, xl_t, yla_t_1, yla_b_1, yl_b] = createiristemplate(eyeimage_filename, pupilx, pupily, pupilr, irisx, irisy, irisr)
% path for writing diagnostic images
global DIAGPATH
DIAGPATH = 'diagnostics';

%normalisation parameters
radial_res = 20;
angular_res = 240;
% with these settings a 9600 bit iris template is
% created

%feature encoding parameters
nscales=1;
minWaveLength=18;
mult=1; % not applicable if using nscales = 1
sigmaOnf=0.5;

%% added

% % eyeimage_filename is a string
% parts = strsplit(eyeimage_filename, '\'); % remove '\'
% 
% % for ND_IRIS 
% filename = parts{11}; % take element 10 (filename) C:\Users\Arif\Documents\MATLAB\masekv3\ND_IRIS\LG4000\images\diagnostic2\04261d1016.tiff
% 
% eyeimage = imread(filename);
eyeimage = imread(eyeimage_filename);

%% 


savefile = [eyeimage_filename,'-houghpara.mat'];
[stat,mess]=fileattrib(savefile);

if stat == 1
    % if this file has been processed before
    % then load the circle parameters and
    % noise information for that file.
    load(savefile);
    
else
    
    % if this file has not been processed before
    % then perform automatic segmentation and
    % save the results to a file
    
    [circleiris, circlepupil, imagewithnoise, xl_t, yla_t, yla_b, yl_b] = segmentiris(eyeimage, pupilx, pupily, pupilr, irisx, irisy, irisr);
%     [circleiris circlepupil imagewithnoise] = segmentiris(eyeimage, eyeimage_filename, pupilx, pupily, pupilr, irisx, irisy, irisr);
%     save(savefile,'circleiris','circlepupil','imagewithnoise');
    
end

% WRITE NOISE IMAGE
%

imagewithnoise2 = uint8(imagewithnoise);
imagewithcircles = uint8(eyeimage);

%get pixel coords for circle around iris
% [x,y] = circlecoords([circleiris(2),circleiris(1)],circleiris(3),size(eyeimage));
[x,y] = circlecoords([cell2mat(irisx),cell2mat(irisy)],cell2mat(irisr),size(eyeimage)); % edited
ind2 = sub2ind(size(eyeimage),double(y),double(x)); 

%get pixel coords for circle around pupil
% [xp,yp] = circlecoords([circlepupil(2),circlepupil(1)],circlepupil(3),size(eyeimage));
% [x,y] = circlecoords([pupily,pupilx],pupilr,size(eyeimage));
[xp,yp] = circlecoords([cell2mat(pupilx),cell2mat(pupily)],cell2mat(pupilr),size(eyeimage)); % edited
ind1 = sub2ind(size(eyeimage),double(yp),double(xp));


% Write noise regions
% imagewithnoise2(ind2) = 255;  % edited
% imagewithnoise2(ind1) = 255;
% Write circles overlayed
imagewithcircles(ind2) = 255;
imagewithcircles(ind1) = 255;
w = cd;
cd(DIAGPATH);
imwrite(imagewithnoise2,[eyeimage_filename,'-noise.jpg'],'jpg');
% imwrite(imagewithcircles,[eyeimage_filename,'-segmented.jpg'],'jpg');
cd(w);

% perform normalisation

% [polar_array noise_array] = normaliseiris(imagewithnoise, circleiris(2),...
%     circleiris(1), circleiris(3), circlepupil(2), circlepupil(1), circlepupil(3),eyeimage_filename, radial_res, angular_res);

[polar_array noise_array] = normaliseiris(imagewithnoise, cell2mat(irisx),...
    cell2mat(irisy), cell2mat(irisr), cell2mat(pupilx), cell2mat(pupily), cell2mat(pupilr),eyeimage_filename, radial_res, angular_res);


% WRITE NORMALISED PATTERN, AND NOISE PATTERN
w = cd;
cd(DIAGPATH);
% imwrite(polar_array,[eyeimage_filename,'-polar.jpg'],'jpg');
% imwrite(noise_array,[eyeimage_filename,'-polarnoise.jpg'],'jpg');
cd(w);

% perform feature encoding
[template mask] = encode(polar_array, noise_array, nscales, minWaveLength, mult, sigmaOnf); 

yla_b_1 = yla_b;
yla_t_1 = yla_t;