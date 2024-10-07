function channelpos = fillmissingprobe(channelpos)
% Fill missing values of a neuropixels probe map

% Extract unique x values
[xUnique,~,chOrdertmp] = unique(channelpos(~isnan(channelpos(:,1)),1),'stable');
xUnique = xUnique(~isnan(xUnique));

chOrder = nan(1,size(channelpos,1));
chOrder(~isnan(channelpos(:,1))) = chOrdertmp;

% Unique y positions per x position
yUnique = arrayfun(@(X) unique(channelpos(channelpos(:,1)==X,2)),xUnique,'Uni',0);

MissingChan = nan(0,2);
ProbablechIdx = nan(0,1);

MissIdx = find(isnan(channelpos(:,1)));
countid = 1;
% Recreate a channelmap based on these values
for xid = 1:numel(xUnique)
    yAllPos = yUnique{xid}; % Current column

    % Most similar column
    SameIdx = find(cellfun(@(X) all(ismember(yUnique{xid},X)) | all(ismember(X,yUnique{xid})),yUnique)); % This column should have the same ypos
    SameIdx(SameIdx == xid) =  [];
    [~,Idx2] = min(abs(numel(yUnique{xid})-cellfun(@numel,yUnique(SameIdx)))); % Most similar column
    yAllPosSim = yUnique{SameIdx(Idx2)};

    if length(unique(diff(yAllPosSim)))>1
        % Also missing a channel, fill in
        yAllPosSim = [min(yAllPosSim):min(diff(yAllPosSim)):max(yAllPosSim)]';
    end

    if any(~ismember(yAllPosSim,yAllPos))
        MissingChan = cat(1,MissingChan,[repmat(xUnique(xid),sum(~ismember(yAllPosSim,yAllPos)),1),yAllPosSim(~ismember(yAllPosSim,yAllPos))]);
        TheseChs = find(chOrder==xid);
        steps = [nan diff(TheseChs)];
        stepsz = nanmin(steps);
        PotentialChans = min(TheseChs)-stepsz*2:stepsz:max(TheseChs);
        PotentialChans(PotentialChans<=0) = [];
        ThisChans = MissIdx(ismember(MissIdx,PotentialChans));
        if isempty(ThisChans)
            ThisChans = MissIdx(countid); % Probably the last channel?
        end
        ProbablechIdx = cat(1,ProbablechIdx,ThisChans);
        countid = countid+1;
    end

end
if any(arrayfun(@(X) sum(ProbablechIdx==X),unique(ProbablechIdx))>1)
    error('Channel map filling in failed')
end
% Compare original and fill in missing channels
channelpos(ProbablechIdx,:) = MissingChan;

return