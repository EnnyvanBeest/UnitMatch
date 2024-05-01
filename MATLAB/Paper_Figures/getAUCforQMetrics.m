function [AUCPerqParam, qParamNames] = getAUCforQMetrics(UMFiles) 
    
    for midx = 1:length(UMFiles)

        fprintf('Reference %s...\n', UMFiles{midx})
        tic
    
        tmpfile = dir(UMFiles{midx});
        if isempty(tmpfile)
            continue
        end

        load(fullfile(tmpfile.folder, tmpfile.name), 'UMparam')

        %% qParams
        noQmetric = 0;
        if ~exist(fullfile(tmpfile.folder, 'qMetricAUCs.mat'))
            try
                QualityMetricsROCs(UMparam.SaveDir)
                close all
            catch ME
                noQmetric = 1;
            end
        end
        if ~noQmetric
            load(fullfile(tmpfile.folder, 'qMetricAUCs.mat'))
    
            if ~exist('qParamNames')
                qParamNames = AUCqParams.qParamNames;
                AUCPerqParam = nan(length(qParamNames),0);
            end
            [takethese,puthere] = ismember(AUCqParams.qParamNames,qParamNames);
            tmpQP = nan(length(qParamNames),1);
            tmpQP(puthere(puthere~=0)) = AUCqParams.AUCMvNM(takethese);
            AUCPerqParam = cat(2,AUCPerqParam,tmpQP);
        end
        toc
    end

    %% AUC distributions for qMetrics
    figure('name','AUC Distr')
    stepsz = 0.05;
    binvec = [0:stepsz:1];
    plotvec = stepsz/2:stepsz:1-stepsz/2;
    for qid = 1:length(qParamNames)
        subplot(ceil(sqrt(length(qParamNames))),round(sqrt(length(qParamNames))),qid)
        s = 1;
        if nanmedian(AUCPerqParam(qid,:))<0.5
            s = -1;
        end
        nums = histcounts(s*AUCPerqParam(qid,:),binvec);
        plot(plotvec,nums,'k-');
        hold on
        line([nanmedian(s*AUCPerqParam(qid,:)) nanmedian(s*AUCPerqParam(qid,:))],get(gca,'ylim'),'color',[1 0 0])
        [h,p(qid)] = ttest(s*AUCPerqParam(qid,:),0.5);

        title([qParamNames{qid} ', p=' num2str(round(p(qid)*100)./100)])
        xlim([0 1])
        makepretty
        offsetAxes
    end