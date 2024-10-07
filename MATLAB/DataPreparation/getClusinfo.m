function clusinfo = getClusinfo(AllKiloSortPaths)

%% Here we're going to actually load in all the sessions requested - only clusinfo to save memory for unitmatch
clusinfo = cell(1, length(AllKiloSortPaths));
addthis = 0;
for subsesid = 1:length(AllKiloSortPaths)
    if isempty(dir(fullfile(AllKiloSortPaths{subsesid}, '*.npy'))) && ~exist(fullfile(AllKiloSortPaths{subsesid}, 'PreparedData.mat'))
        continue
    end

    disp(['Loading clusinfo for ', AllKiloSortPaths{subsesid}])
    tmp = matfile(fullfile(AllKiloSortPaths{subsesid}, 'PreparedData.mat'));
    clusinfo{subsesid} = tmp.clusinfo;

    % Replace recsesid with subsesid
    clusinfo{subsesid}.RecSesID = clusinfo{subsesid}.RecSesID + addthis;
    addthis = max(clusinfo{subsesid}.RecSesID);
end

% Add all cluster information in one 'cluster' struct - can be used for further analysis
% Apparently sometimes clusinfo does not have the same fields (KS version?)
allfieldNames = cellfun(@fieldnames,clusinfo,'Uni',0);
[unqFields,~,unqID] = unique(vertcat(allfieldNames{:}));
% find common fields
match_ID = find(accumarray(unqID(:),1) == numel(allfieldNames));
commonFields = unqFields(match_ID);
% remove unshared fields
for sesid = 1:length(clusinfo)
    thesefields = fieldnames(clusinfo{sesid});
    clusinfo{sesid} = rmfield(clusinfo{sesid},thesefields(~ismember(thesefields,commonFields))); % REmove fields
end


% Continue adding them together
clusinfo = [clusinfo{:}];
clusinfoNew = struct;
for fieldid = 1:length(commonFields)
    try
        eval(['clusinfoNew.', commonFields{fieldid}, '= cat(1,clusinfo(:).', commonFields{fieldid} ');'])
    catch ME
        try
            eval(['clusinfoNew.', commonFields{fieldid}, '= cat(2,clusinfo(:).', commonFields{fieldid}, ');'])
        catch ME
        end
    end
end
clusinfo = clusinfoNew;
clear clusinfoNew

return
