function [x1, y2] = getpeaks(peaks,Inoise,k1,k2)

% figure, imshow(Inoise);
% 
% hold on;
for peak = peaks
    [x, y] = circlepoints(peak(3));
    x1 = x+peak(1);
    y2 = y+peak(2);
%     if (y2 >= k1 && y2 <= k2)
%     plot(x+peak(1), y+peak(2), 'g-');
%         plot(x1, y2, 'g-');
%     end
end

% hold off

end