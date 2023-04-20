function ExtractFields(StructNames)

% Extracts the parameters from struct(s) StructNames (cell with the
% struct names (struct, not characters!))
for sid = 1:length(StructNames)
    tmp = StructNames{sid};
    FN = fieldnames(tmp);
    for fid=1:length(FN)
        eval([FN{fid} '= tmp.' FN{fid} ';'])
        assignin('caller',FN{fid},eval(FN{fid}))
    end
end

if exist('clu','var')
    % Asign spike shank
    spikeShank = nan(length(clu),1);
    ShankOpt = unique(Shank);
    for shid = 1:length(ShankOpt)
        spikeShank(ismember(clu,cluster_id(Shank==ShankOpt(shid)))&ismember(RecSes,RecSesID(Shank==ShankOpt(shid)))) = ShankOpt(shid);
    end
    assignin('caller','spikeShank',spikeShank)
end
