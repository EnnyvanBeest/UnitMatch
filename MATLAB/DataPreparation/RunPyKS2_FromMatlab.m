%% Automated
% Load all data
% Find available datasets (always using dates as folders)
clear DateOpt
DateOpt = arrayfun(@(X) dir(fullfile(DataDir{DataDir2Use(X)},MiceOpt{X},'*-*')),1:length(MiceOpt),'UniformOutput',0);
DateOpt = cellfun(@(X) X([X.isdir]),DateOpt,'UniformOutput',0);
DateOpt = cellfun(@(X) {X.name},DateOpt,'UniformOutput',0);

maxses = inf; %For debugging we might want to do only maxses nr at a time
for midx = 1:length(MiceOpt)
    myKsDir = fullfile(KilosortDir,MiceOpt{midx});
    subksdirs = dir(fullfile(myKsDir,'Probe*')); %This changed because now I suddenly had 2 probes per recording

    if strcmp(RecordingType{midx},'Chronic') && PipelineParams.RunPyKSChronicStitched % These are my chronic mice and we want pyks to do the matching, one dataset per mouse
        countid=1;
        % For every date a different dataset
        Dates4Mouse = DateOpt{midx};
        IMROTable = cell(1,length(Dates4Mouse)); %Save out IMRO table used
        FullPath = cell(1,length(Dates4Mouse)); % Save out full path name for later use
        ProbeID = cell(1,length(Dates4Mouse)); %Probe ID in case there's multiple probes

        TakeIMRO = {};

        for didx = 1:length(Dates4Mouse)
            % Within folders, look for 'RF mapping sessions'
            thisdate = Dates4Mouse{didx};

            %% Loading data from kilosort/phy easily
            myKsDir = fullfile(KilosortDir,MiceOpt{midx},thisdate);
            tmpephysdir = dir(fullfile(DataDir{DataDir2Use(midx)},MiceOpt{midx},thisdate,'ephys',['*' MiceOpt{midx} '*']));
            if exist('IgnoreTheseFiles','var')
                for id = 1:length(IgnoreTheseFiles)
                    % Check if there's an ephys folder, if so run pyks2
                    tmpephysdir(find(cell2mat(cellfun(@(X) any(strfind(X,IgnoreTheseFiles{id})),{tmpephysdir(:).name},'UniformOutput',0))))=[];
                end
            end
            if isempty(tmpephysdir)
                continue
            end

            % Copy file and then run pykilosort on it
            for id = 1:length(tmpephysdir)
                tmpfile = dir(fullfile(tmpephysdir(id).folder,tmpephysdir(id).name,'**','*ap*bin'));

               
                if isempty(tmpfile)
                    continue
                end
                for sesid = 1:length(tmpfile)
                    ProbeName = strsplit(tmpfile(sesid).name,'imec');
                    ProbeName = ['Probe' num2str(ProbeName{2}(1))];

                    Sesinfo = strsplit(tmpephysdir(id).name,'_g0');
                    Sesinfo = strsplit(Sesinfo{1},MiceOpt{midx});
                    Sesinfo = Sesinfo{2};
                    if isempty(Sesinfo)
                        Sesinfo = '';
                    else
                        Sesinfo = Sesinfo(end);
                    end
                    % Read meta file
                    meta = ReadMeta2(fullfile(tmpfile(sesid).folder));
                    IMROTable{countid}=meta.imroTbl;
                    FullPath{countid}= fullfile(tmpfile(sesid).folder,tmpfile(sesid).name);
                    ProbeID{countid} = ProbeName;
                    countid=countid+1;

                end
            end
        end
        IMROTable(cellfun(@isempty,IMROTable))=[];
        IMROTableOpt = unique(IMROTable); %Unique IMRO tables
        ProbeID(cellfun(@isempty,ProbeID))=[];
        ProbeIDOpt = unique(ProbeID);
        FullPath(cellfun(@isempty,FullPath)) = [];

        for pid = 1:length(ProbeIDOpt)
            for id = 1:length(IMROTableOpt)
                % Define save path
                SaveDirTmp = fullfile(KilosortDir, MiceOpt{midx},...
                    'Chronic',ProbeIDOpt{pid},['IMRO_' num2str(id)]);


                theseses = find(ismember(IMROTable,IMROTableOpt{id})&ismember(ProbeID,ProbeIDOpt{pid})); % All sessions with this IMRO table
                ThesePaths = FullPath(theseses);


                disp('Check whether this chronic data set already exists')
                doneflag = 0;
                if ~exist(SaveDirTmp)
                    mkdir(SaveDirTmp)
                else
                    FileList = dir(fullfile(fullfile(KilosortDir, MiceOpt{midx},...
                        'Chronic',ProbeIDOpt{pid}),'*','SessionsIncluded.mat'));
                    if ~isempty(FileList)
                        for fid = 1:length(FileList)
                            alreadydoneses = load(fullfile(FileList(fid).folder,FileList(fid).name));
                            if ~any(~ismember(alreadydoneses.ThesePaths,ThesePaths))
                                TakeIMRO = {TakeIMRO{:} FileList(fid).folder}
                                disp('Already done... skip')
                                doneflag=1;
                                %                             continue
                            end
                        end
                    end
                end
                if doneflag
                    continue
                end

                if ~exist(fullfile(tmpdatafolder,['IMRO_' num2str(id)]))
                    mkdir(fullfile(tmpdatafolder,['IMRO_' num2str(id)]))
                end
                if length(ThesePaths)>maxses
                    theseses = theseses(1:maxses)
                    ThesePaths = ThesePaths(1:maxses);
                end
                % loop over included paths and save data locally
                dat_path = ['['];
                for sesid = 1:length(ThesePaths)
                    tmpfile = dir(ThesePaths{sesid});
                    if ~exist(fullfile(tmpdatafolder,['IMRO_' num2str(id)],tmpfile(1).name)) & ~doneflag
                        disp('Copying data to local folder...')
                        copyfile(fullfile(tmpfile(1).folder,tmpfile(1).name),fullfile(tmpdatafolder,['IMRO_' num2str(id)],tmpfile(1).name))
                        metafile = dir(fullfile(tmpfile(1).folder,'*ap.meta'));
                        copyfile(fullfile(metafile.folder,metafile.name),fullfile(tmpdatafolder,['IMRO_' num2str(id)],metafile.name))
                        chfile = dir(fullfile(tmpfile(1).folder,'*ap.ch'));
                        copyfile(fullfile(chfile.folder,chfile.name),fullfile(tmpdatafolder,['IMRO_' num2str(id)],chfile.name))
                    end
                    dat_path = [dat_path '"' strrep(fullfile(tmpfile(1).folder,tmpfile(1).name),'\','/') '", '];
                end
                dat_path(end-1)=[];
                dat_path(end) = ']'
                % PyKS2
                if ~doneflag
                    RunFromFolder = fullfile(tmpdatafolder,['IMRO_' num2str(id)]);                   
                    success = pyrunfile("RunPyKS2_FileInput.py","success",ThisFile = strrep(fullfile(RunFromFolder),'\','/'))
                    clear success

                    % now copy the output
                    disp('Copying output from temporary directory to KS folder')
                    try
                        copyfile(fullfile(RunFromFolder,['IMRO_' num2str(id)],'output','*'),fullfile(SaveDirTmp))
                    catch
                        movefile(fullfile(RunFromFolder,'output','*'),fullfile(SaveDirTmp))
                    end

                    % Remove temporary files
                    disp('Delete temporary folder')
                    try
                        subfolds = dir(fullfile(tmpdatafolder,['IMRO_' num2str(id), '\**']))
                        arrayfun(@(X) delete(fullfile(X.folder,X.name)),subfolds,'UniformOutput',0)
                        subfolds = dir(fullfile(tmpdatafolder,['IMRO_' num2str(id), '\**']))
                        subfolds = unique({subfolds(:).folder});
                        subfolds = fliplr(sort(subfolds));
                        cellfun(@(X) rmdir(X),subfolds,'UniformOutput',0)
                    catch ME
                        disp(ME)
                    end

                    % Save which files and the IMRO table
                    save(fullfile(SaveDirTmp,'SessionsIncluded.mat'),'ThesePaths')

                    channelpostmpconv = ChannelIMROConversion(fileparts(FullPath{theseses(1)}),1); % For overview of which channels were recorded
                    drawnow
                    saveas(gcf,fullfile(SaveDirTmp,'IMRO.fig'))
                    saveas(gcf,fullfile(SaveDirTmp,'IMRO.bmp'))

                end
   
                % change dat_path in params file
                paramsfile = dir(fullfile(SaveDirTmp,'params.py'));
                fid =fopen(fullfile(paramsfile.folder,paramsfile.name));
                C=textscan(fid,'%s','delimiter','\n');
                fclose(fid);
                idx = find(cell2mat(cellfun(@(X) any(strfind(X,'dat_path')),C{1},'UniformOutput',0)));
                
                tmp = regexp(C{1}{idx},'='); %Find =
                C{1}{idx}(tmp+1:end) = []; %Remove current paths
                C{1}{idx} = strcat(C{1}{idx}, dat_path)
               
                % Rename old params
                movefile(fullfile(SaveDirTmp,'params.py'),fullfile(SaveDirTmp,'paramsOri.py'))
                % print new file
                fName = fullfile(SaveDirTmp,'params.py');
                fid = fopen(fName,'w');            % Open the file
                for k=1:numel(C{1,1})
                    fprintf(fid,'%s\r\n',C{1,1}{k,1});
                end
                fclose(fid);

            end

        end
    else
        % For every date a different dataset
        Dates4Mouse = DateOpt{midx};
        for didx = 1:length(Dates4Mouse)
            % Within folders, look for 'RF mapping sessions'
            thisdate = Dates4Mouse{didx};
            %% Loading data from kilosort/phy easily
            myKsDir = fullfile(KilosortDir,MiceOpt{midx},thisdate);

            % Check if there's an ephys folder, if so run pyks2
            tmpephysdir = dir(fullfile(DataDir{DataDir2Use(midx)},MiceOpt{midx},thisdate,'ephys',['*' MiceOpt{midx} '*']));
            if isempty(tmpephysdir)
                continue
            end
            if exist('IgnoreTheseFiles','var')
                for id = 1:length(IgnoreTheseFiles)
                    % Check if there's an ephys folder, if so run pyks2
                    tmpephysdir(find(cell2mat(cellfun(@(X) any(strfind(X,IgnoreTheseFiles{id})),{tmpephysdir(:).name},'UniformOutput',0))))=[];
                end
            end

            if isempty(tmpephysdir)
                continue
            end

            % Copy file and then run pykilosort on it
            for id = 1:length(tmpephysdir)
                tmpfile = dir(fullfile(tmpephysdir(id).folder,tmpephysdir(id).name,'**','*ap*bin'));


                for sesid = 1:length(tmpfile)

                    if isempty(tmpfile) %|| any(strfind(tmpfile(sesid).folder,'NatImages')) ||  any(strfind(tmpfile(sesid).folder,'Spontaneous')) % Not my session
                        continue
                    end
                    ProbeName = strsplit(tmpfile(sesid).name,'imec');
                    ProbeName = ['Probe' num2str(ProbeName{2}(1))];

                    Sesinfo = strsplit(tmpephysdir(id).name,'_g0');
                    Sesinfo = strsplit(Sesinfo{1},MiceOpt{midx});
                    Sesinfo = Sesinfo{2};
                    if isempty(Sesinfo)
                        Sesinfo = strsplit(tmpephysdir(id).name,'_g0');
                        Sesinfo = strsplit(Sesinfo{2},'_');
                        Sesinfo = Sesinfo{2};
                        if isempty(Sesinfo)

                        Sesinfo = '';
                        end
                    else
                        Sesinfo = Sesinfo(end);
                    end

                    %                     myKsDir = fullfile(myKsDir,['Probe' num2str(probeid-1)])
                    subksdirs = dir(fullfile(myKsDir,ProbeName,Sesinfo)); %This changed because now I suddenly had 2 probes per recording

                    if ~isempty(subksdirs)
                        continue
                    end

                    if ~isempty(dir(fullfile(myKsDir,ProbeName,Sesinfo,'spike*')))
                        disp([fullfile(myKsDir,ProbeName,Sesinfo) ' already done, skip...'])
                        continue
                    end
                    metafile = dir(fullfile(tmpfile(sesid).folder,'*ap.meta'));
                    try
                        chfile = dir(fullfile(tmpfile(sesid).folder,'*ap.ch'));
                        if ~isempty(chfile)
                            Compressed=1;
                        else
                            Compressed=0;
                        end

                    catch ME
                        Compressed=0;
                    end

                    if ~exist(fullfile(tmpdatafolder,tmpfile(sesid).name))
                        disp('Copying data to local folder...')
                        copyfile(fullfile(tmpfile(sesid).folder,tmpfile(sesid).name),fullfile(tmpdatafolder,tmpfile(sesid).name))
                        copyfile(fullfile(metafile.folder,metafile.name),fullfile(tmpdatafolder,metafile.name))
                        if Compressed
                        copyfile(fullfile(chfile.folder,chfile.name),fullfile(tmpdatafolder,chfile.name))
                        end
                    end
                    % PyKS2
                    try
                        success = pyrunfile("RunPyKS2_FromMatlab.py","success",bin_file = strrep(fullfile(tmpdatafolder,tmpfile(sesid).name),'\','/'))
                        clear success
                    catch ME
                        disp(ME)
                        disp([fullfile(tmpdatafolder,tmpfile(sesid).name) ' not successfull... skip'])
                        continue
                    end

                    % now copy the output
                    disp('Copying output from temporary directory to KS folder')
                    copyfile(fullfile(tmpdatafolder,'output'),fullfile(myKsDir,ProbeName,Sesinfo))

                    % Remove temporary files
                    disp('Delete files from temporary folder')
                    try
                        delete(fullfile(tmpdatafolder,metafile.name))
                        delete(fullfile(tmpdatafolder,chfile.name))
                    catch ME
                        disp(ME)
                    end
                    try
                        delete(fullfile(tmpdatafolder,tmpfile(sesid).name))
                    catch ME
                        disp(ME)
                    end
                    try
                        delete(fullfile(tmpdatafolder,'output','*'))
                        rmdir(fullfile(tmpdatafolder,'output'))
                    catch ME
                        disp(ME)
                    end

                    % change dat_path in params file
                    paramsfile = dir(fullfile(myKsDir,ProbeName,Sesinfo,'params.py'));
                    fid =fopen(fullfile(paramsfile.folder,paramsfile.name));
                    C=textscan(fid,'%s','delimiter','\n');
                    fclose(fid);
                    idx = find(cell2mat(cellfun(@(X) any(strfind(X,'dat_path')),C{1},'UniformOutput',0)));

                    tmp = regexp(C{1}{idx},'='); %Find =
                    C{1}{idx}(tmp+1:end) = []; %Remove current paths
                    C{1}{idx} = strcat(C{1}{idx}, ['r"' strrep(fullfile(tmpfile(sesid).folder,tmpfile(sesid).name),'\','/') '"'])

                    % Rename old params
                    movefile(fullfile(myKsDir,ProbeName,Sesinfo,'params.py'),fullfile(myKsDir,ProbeName,Sesinfo,'paramsOri.py'))
                    % print new file
                    fName = fullfile(myKsDir,ProbeName,Sesinfo,'params.py');
                    fid = fopen(fName,'w');            % Open the file
                    for k=1:numel(C{1,1})
                        fprintf(fid,'%s\r\n',C{1,1}{k,1});
                    end
                    fclose(fid);

                end
            end


        end
    end
end