% created on 261218

function [x1, y2] = getpeaks3(peaks,I,eyeimage_filename)

figure, imshow(I);

hold on;

for peak = peaks
    [x, y] = circlepoints(peak(3));
    x1 = x+peak(1);
    y2 = y+peak(2);
%     if (y2 >= k1 && y2 <= k2)
%     plot(x+peak(1), y+peak(2), 'g-');
        plot(x1, y2, 'g-');
%     end

s=getframe;
    
end

hold off
% 
% w = cd;
% cd(DIAGPATH);

imwrite(s.cdata,[eyeimage_filename,'.tiff-homomorphic-lens_segmented.jpg'],'jpg');

close all;

% cd(w);
end