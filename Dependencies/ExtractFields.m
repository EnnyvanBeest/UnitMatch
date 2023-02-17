function ExtractFields(StructNames)

% Extracts the parameters from struct(s) StructNames (cell with the
% struct names (struct, not characters!))
for sid = 1:length(StructNames)
    tmp = StructNames{sid};
    FN = fieldnames(tmp);
    for fid=1:length(FN)
        eval([FN{fid} '= tmp.' FN{fid} ';'])
        assignin('base',FN{fid},eval(FN{fid}))
    end
end

% Asign spike shank
spikeShank = nan(length(clu),1);
ShankOpt = unique(Shank);
for shid = 1:length(ShankOpt)
    spikeShank(ismember(clu,cluster_id(Shank==ShankOpt(shid)))&ismember(RecSes,RecSesID(Shank==ShankOpt(shid)))) = ShankOpt(shid);
end
assignin('base','spikeShank',spikeShank)

