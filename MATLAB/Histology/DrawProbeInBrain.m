
%% Plot probes
disp(fullfile(SaveDir,MiceOpt{midx},thisdate,thisprobe,[num2str(SN) '_HistoEphysAlignment.mat']))
Depth2AreaPer = load(fullfile(SaveDir,MiceOpt{midx},thisdate,thisprobe,[num2str(SN) '_HistoEphysAlignment.mat']));
Depth2AreaPer=Depth2AreaPer.Depth2Area;

probePoints = cell2mat(Depth2AreaPer.Coordinates)./10; %Because we use 10 micron steps
probePoints(any(isnan(probePoints),2),:)=[];
% Orientation is ASR, has to go ASL
probePoints(:,3) = size(av,3)-probePoints(:,3);

% figure(ProbeFigure)
% hold on
% plot3(probePoints(:,3),probePoints(:,1),probePoints(:,2),'.')
% 
% xlabel('LR')
% zlabel('SI')
% ylabel('AP')
% set(gca,'zdir','reverse')

if any(probePoints)<0
    keyboard
end
% get the probe points for the currently analyzed probe
if strcmp(plane,'coronal')
    curr_probePoints = probePoints(:, [1 2 3]);
elseif strcmp(plane,'sagittal')
    curr_probePoints = probePoints(:, [3 2 1]);
elseif strcmp(plane,'transverse')
    curr_probePoints = probePoints(:, [3 1 2]);
end


% get line of best fit through points
% m is the mean value of each dimension; p is the eigenvector for largest eigenvalue
[m,p,s] = best_fit_line(curr_probePoints(:,1), curr_probePoints(:,2), curr_probePoints(:,3));
if isnan(m(1))
    disp(['no points found for probe '])
else

% ensure proper orientation: want 0 at the top of the brain and positive distance goes down into the brain
if p(2)<0
    p = -p;
end

% determine "origin" at top of brain -- step upwards along tract direction until tip of brain / past cortex
out_of_brain = false;
while ~(ann==1 && out_of_brain) % && distance_stepped > .5*active_probe_length)
    m = m-p; % step 10um, backwards up the track
    ann = av(round(m(1)),round(m(2)),round(m(3))); %until hitting the top
    if strcmp(st_allen.safe_name(ann), 'root')
        % make sure this isn't just a 'root' area within the brain
        m_further_up = m - p*20; % is there more brain 200 microns up along the track?
        ann_further_up = av(round(max(1,m_further_up(1))),round(max(1,m_further_up(2))),round(max(1,m_further_up(3))));
        if strcmp(st_allen.safe_name(ann_further_up), 'root')
            out_of_brain = true;
        end
    end
end

% focus on wireframe plot
figure(fwireframe);

% plot line the length of the entire probe in reference space
dist = norm(m-curr_probePoints(end,:));
plot3(m(1)+p(1)*[1 dist], m(3)+p(3)*[1 dist], m(2)+p(2)*[1 dist], ...
    'Color', ProbeColors(midx,:), 'LineWidth', 2, 'LineStyle','-');

% plot probe points
hp = plot3(curr_probePoints(:,1), curr_probePoints(:,3), curr_probePoints(:,2), '.','linewidth',2, 'color',[ProbeColors(midx,:) .2],'markers',10);

% plot brain entry point
plot3(m(1), m(3), m(2), 'r*','linewidth',1)


end
