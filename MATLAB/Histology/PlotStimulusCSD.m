function fig = PlotStimulusCSD(C, titleStr)
% fig = PlotStimulusCSD(C, titleStr)
% Stimulus-evoked LFP and CSD (StimulusCSD output) per shank: LFP (depth x time), CSD (sinks blue,
% sources red, in baseline SDs per depth) with the response onset and the early-sink window, and
% the early sink per depth.
if nargin < 2, titleStr = ''; end
fig = figure('Name','Stimulus CSD','Position',[50 50 min(1800,max(800,650*numel(C))) 750]);
x = linspace(-1,1,256)'; cmapCSD = [min(1,1+x) 1-abs(x) min(1,1-x)]; % blue (sink) - white - red (source)
nS = numel(C);
for k = 1:nS
    c = C(k);
    ax = subplot(1,3*nS,3*k-2);
    imagesc(ax,c.t*1000,c.Depth,c.LFP); set(ax,'YDir','normal'); colormap(ax,parula)
    xline(ax,0,'k:'); title(ax,{sprintf('shank %d',c.Shank),'LFP'}); xlabel(ax,'ms'); ylabel(ax,'depth on probe (um)')
    ax = subplot(1,3*nS,3*k-1);
    im = imagesc(ax,c.t*1000,c.Depth,c.CSD); set(ax,'YDir','normal'); im.AlphaData = ~isnan(c.CSD);
    lim = max(3,prctile(abs(c.CSD(:)),99)); clim(ax,[-lim lim]); colormap(ax,cmapCSD)
    xline(ax,0,'k:'); xline(ax,1000*c.SinkWin,'k--'); yline(ax,c.SinkDepth,'g-','LineWidth',1)
    title(ax,{sprintf('CSD, %d trials',c.nTrials),sprintf('onset %.0f ms',c.OnsetMs)}); xlabel(ax,'ms')
    ax = subplot(1,3*nS,3*k);
    plot(ax,c.EarlySink,c.Depth,'k'); grid(ax,'on'); ylim(ax,c.Depth([1 end]))
    xlabel(ax,'early sink (SD)'); title(ax,{'early sink',sprintf('%.1f SD at %.0f um',c.SinkZ,c.SinkDepth)})
end
if ~isempty(titleStr), sgtitle(titleStr,'Interpreter','none'); end
end
