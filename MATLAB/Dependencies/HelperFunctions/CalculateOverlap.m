function overlap = CalculateOverlap(hist1,hist2,drawnow)
% Calculates the overlap between histogram 1 and histogram 2
if length(hist1)~= length(hist2)
    error('histograms are not equal in length. Are these counts for the same range of x?')
end
if nargin<3
    drawnow=0;
end
d = hist2-hist1; % Normalized sum

if drawnow
    figure; plot(hist1); hold on; plot(hist2); plot(d)
end

overlap = 1-sum(abs(d))./(sum(hist1)+sum(hist2));
if overlap<0
    keyboard
end

end