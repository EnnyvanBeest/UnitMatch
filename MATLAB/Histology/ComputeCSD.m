function C = ComputeCSD(Avg, Shank, Depth, t, Opt)
% C = ComputeCSD(Avg, Shank, Depth, t, Opt)
% Current source density per shank from an event-averaged LFP.
%  Avg: nChan x nT averaged LFP (raw units); Shank (0-based) and Depth (um) per channel; t: nT times (s)
%  Per shank: channels at the same depth (the two columns) are combined by the median, the profile is
%  interpolated on a regular grid (Opt.GridUm), smoothed along depth (Gaussian, sigma Opt.SmoothUm) and
%  CSD = -(phi(z+h) + phi(z-h) - 2 phi(z))/h^2 with h = Opt.DiffUm (arbitrary units; sinks negative).
%  Opt.BaseWin: window (s) for the noise estimate (std of the CSD there, median over depths).
% C: struct array per shank: Shank, Depth (grid, um), t, LFP, CSD (nDepth x nT), Noise
if nargin < 5, Opt = struct; end
Def = struct('GridUm',10,'SmoothUm',30,'DiffUm',40,'BaseWin',[-Inf 0]);
for fn = fieldnames(Def)', if ~isfield(Opt,fn{1}), Opt.(fn{1}) = Def.(fn{1}); end, end
C = struct('Shank',{},'Depth',{},'t',{},'LFP',{},'CSD',{},'Noise',{});
t = t(:)';
for sh = unique(Shank(:))'
    ch = find(Shank == sh);
    if numel(ch) < 10, continue; end
    [ud,~,g] = unique(round(Depth(ch)/5)*5); % same depth within 5 um
    if numel(ud) < 10, continue; end
    P = zeros(numel(ud),numel(t));
    for k = 1:numel(ud), P(k,:) = median(Avg(ch(g == k),:),1); end
    grid = (ceil(ud(1)/Opt.GridUm)*Opt.GridUm : Opt.GridUm : floor(ud(end)/Opt.GridUm)*Opt.GridUm)';
    G = interp1(ud,P,grid,'linear');
    G = smoothdata(G,1,'gaussian',max(3,round(5*Opt.SmoothUm/Opt.GridUm)));
    h = max(1,round(Opt.DiffUm/Opt.GridUm));
    csd = nan(size(G));
    csd(1+h:end-h,:) = -(G(1+2*h:end,:) + G(1:end-2*h,:) - 2*G(1+h:end-h,:))/(h*Opt.GridUm/100)^2;
    % the smoothing is one-sided near the ends of the probe: no CSD there
    edge = h + round(2*Opt.SmoothUm/Opt.GridUm);
    csd([1:min(edge,end) max(1,end-edge+1):end],:) = NaN;
    base = t >= Opt.BaseWin(1) & t < Opt.BaseWin(2);
    noise = median(std(csd(:,base),0,2,'omitnan'),'omitnan');
    C(end+1) = struct('Shank',sh,'Depth',grid,'t',t,'LFP',G,'CSD',csd,'Noise',noise); %#ok<AGROW>
end
end
