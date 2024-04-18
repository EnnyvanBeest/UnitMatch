function channelpos = fillmissingprobe(channelpos)
% Fill missing values of a neuropixels probe map

% Extract unique values
xUnique = unique(channelpos(:,1),'stable');
yUnique = unique(channelpos(:,2),'stable');

xUnique = xUnique(~isnan(xUnique));
yUnique = yUnique(~isnan(yUnique));

% Steps in y
ySteps = unique(diff(yUnique));

% identify nans
rowId = find(isnan(channelpos(:,1)));
for id = 1:numel(rowId)
    % Take a few samples around the missing channel
    TakeRows = rowId-length(xUnique)/2:rowId+length(xUnique)/2; 
    TakeRows(TakeRows<1|TakeRows>size(channelpos,1)) = [];

    % Extract piece of channel
    channelpostmp = channelpos(TakeRows,:);

    % pattern completion
    pattern = arrayfun(@(X) find(channelpostmp(:,1)==X,1,'first'),xUnique,'Uni',0);
    MissingVal = xUnique(cellfun(@isempty,pattern));
    channelpos(rowId,1) = MissingVal;

    % Now for the second dimension
    pattern = arrayfun(@(X) find(channelpostmp(:,2)==X,1,'first'),yUnique,'Uni',0);
    usedyVals = yUnique(~cellfun(@isempty,pattern));
    [~,id2] = ismember(channelpostmp(:,2),usedyVals);
    if mod(rowId,2) == 0 % even; repeat old number
        channelpos(rowId,2) = usedyVals(id2(find(id2==0)-1));
    else % start new number
        channelpos(rowId,2) = usedyVals(id2(find(id2==0)-1)) + ySteps;
    end
end






return