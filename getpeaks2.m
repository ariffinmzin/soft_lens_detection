function [x1, y2] = getpeaks2(peaks,Inoise,k1,k2,xx,eyeimage_filename)

global DIAGPATH;

% k1
% k2
aa = [k1 , k2];
l = length(xx);

figure, imshow(Inoise);

hold on;



% for peak = peaks % asal
for j = 1:l
%     peak
    [x, y] = circlepoints(peaks(3));
    x1(j) = x(j)+peaks(1);
    y2(j) = y(j)+peaks(2);
    if (y2(j) >= aa(1) && y2(j) <= aa(2))
        yy2(j) = y2(j); 
        yy2(yy2 == 0) = NaN;  
        xx1(j) = x1(j);
        xx1(xx1 == 0) = NaN;
        plot(xx1, yy2, 'g-');

    end
    
    s=getframe;

end

%%%

hold off
% 
w = cd;
cd(DIAGPATH);

imwrite(s.cdata,[eyeimage_filename,'-lens_segmented.jpg'],'jpg');

close all;

cd(w);

%%%











% f = figure('visible', 'off');
% hold on;
% 
% for j = 1:l
% %     peak
%     [x, y] = circlepoints(peaks(3));
%     x1(j) = x(j)+peaks(1);
%     y2(j) = y(j)+peaks(2);
%     if (y2(j) >= aa(1) && y2(j) <= aa(2))
%         yy2(j) = y2(j); 
%         yy2(yy2 == 0) = NaN;  
%         xx1(j) = x1(j);
%         xx1(xx1 == 0) = NaN;
%         plot(xx1, yy2, 'g-');
% 
%     end
%     
%     s=getframe;
% 
% end
% 
% hold off;
% close(f)
% 
% imwrite(s.cdata,[eyeimage_filename,'-lens_segmented.jpg'],'jpg');

end