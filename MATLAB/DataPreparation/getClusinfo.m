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
fields = fieldnames(clusinfo{1});
IncludeField = true(1,length(fields));
for sesid=2:length(clusinfo)
    IncludeField(~ismember(fields,fieldnames(clusinfo{sesid}))) = false;
end
for sesid=1:length(clusinfo)
    thesefields = fieldnames(clusinfo{sesid});
    clusinfo{sesid} = rmfield(clusinfo{sesid},thesefields(ismember(thesefields,fields(find(~IncludeField))))); % REmove fields
end

% Continue adding them together
clusinfo = [clusinfo{:}];
clusinfoNew = struct;
fields = {fields{find(IncludeField)}}';
for fieldid = 1:length(fields)
    try
        eval(['clusinfoNew.', fields{fieldid}, '= cat(1,clusinfo(:).', fields{fieldid} ');'])
    catch ME
        try
            eval(['clusinfoNew.', fields{fieldid}, '= cat(2,clusinfo(:).', fields{fieldid}, ');'])
        catch ME
            keyboard
        end
    end
end
clusinfo = clusinfoNew;
clear clusinfoNew

return
