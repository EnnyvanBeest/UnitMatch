% Mixture of normal distributions
function [MSE yPred] = MixtureOfNorms(Coefficients,x,y)
% inputs:
% Coefficients with [mu std frac, mu std frac], etc. per distribution
% x: x-values
nNorm = length(Coefficients)./3;
yPred = zeros(1,length(x));
for normid = 1:nNorm
    tmpnorm = normpdf(x,Coefficients((normid-1)*3+1),Coefficients((normid-1)*3+2));
    tmpnorm = (tmpnorm-min(tmpnorm(:)))./(max(tmpnorm(:))-min(tmpnorm(:)));
    yPred = yPred + Coefficients(normid*3).*tmpnorm;
end
MSE = nanmean((y-yPred).^2);

if 1
    tmpfig = figure;
    plot(x,y,'k.')
    hold on
    plot(x,yPred,'r-')
    title(num2str(MSE))
    pause(1)
    close(tmpfig)
end
return