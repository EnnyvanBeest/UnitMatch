
%% Tracked neurons
% Liberal
DateOpt = cellfun(@(X) strsplit(X,'\'),UMparam.KSDir,'Uni',0);
DateOpt = cellfun(@(X) datetime(X{5}),DateOpt,'Uni',0);
DateInDeltaDays = cell2mat(cellfun(@(X) days(X-DateOpt{1}),DateOpt,'Uni',0));

DayIdx = UniqueIDConversion.recsesAll(logical(UniqueIDConversion.GoodID));
DayDelta = DateInDeltaDays(DayIdx);
AllUIDs = (UniqueIDConversion.UniqueID(logical(UniqueIDConversion.GoodID)));
AllUIDsCons = (UniqueIDConversion.UniqueIDConservative(logical(UniqueIDConversion.GoodID)));

DayOpt = unique(DateInDeltaDays);
PercStableLiberal = nan(3,0); % DifferenceInDays X PercStable x drift
PercStableConservative = nan(3,0); % DifferenceInDays X PercStable x drift

for did1 = 1:length(DayOpt)
    NeuronsD1 = AllUIDsCons(DayDelta==DayOpt(did1));
    NeuronsD1Cons =  AllUIDs(DayDelta==DayOpt(did1));
    for did2 = did1:length(DayOpt)

        % What was the Z drift?
        tmpdrift = abs(nanmean(UMparam.drift(DateInDeltaDays==DayOpt(did2),3,1)) - nanmean(UMparam.drift(DateInDeltaDays==DayOpt(did1),3,1)))

        NeuronsD2 = AllUIDs(DayDelta==DayOpt(did2));
        PercStableLiberal = cat(2,PercStableLiberal,[DayOpt(did2)-DayOpt(did1),sum(ismember(NeuronsD1,NeuronsD2))/min([numel(NeuronsD1),numel(NeuronsD2)]),tmpdrift]');

        NeuronsD2Cons = AllUIDsCons(DayDelta==DayOpt(did2));
        PercStableConservative = cat(2,PercStableConservative,[DayOpt(did2)-DayOpt(did1),sum(ismember(NeuronsD1Cons,NeuronsD2Cons))/min([numel(NeuronsD1Cons),numel(NeuronsD2Cons)]),tmpdrift]');

    end
end
figure('name','Tracking Across Days')
subplot(2,2,1)
scatter(PercStableLiberal(1,:),PercStableLiberal(2,:),20,[0 0 0],'filled')
hold on
scatter(PercStableConservative(1,:),PercStableConservative(2,:),20,[1 0 0],'filled')
ylabel('Fraction tracked neurons')
xlabel('\Delta days')
title('Tracking across days')
makepretty
offsetAxes

subplot(2,2,2)
hold on
Edges = min(PercStableLiberal(1,:)):10:max(PercStableLiberal(1,:));
Vals =  min(PercStableLiberal(1,:))+5:10:max(PercStableLiberal(1,:))-5;
ValPerGroupMean = cell2mat(arrayfun(@(X) nanmean(PercStableLiberal(2,(PercStableLiberal(1,:)>=Edges(X)&PercStableLiberal(1,:)<=Edges(X+1)))),1:length(Edges)-1,'Uni',0));
ValPerGroupstd = cell2mat(arrayfun(@(X) nanstd(PercStableLiberal(2,(PercStableLiberal(1,:)>=Edges(X)&PercStableLiberal(1,:)<=Edges(X+1)))),1:length(Edges)-1,'Uni',0));
ValPerGroupMean2 = cell2mat(arrayfun(@(X) nanmean(PercStableConservative(2,(PercStableConservative(1,:)>=Edges(X)&PercStableConservative(1,:)<=Edges(X+1)))),1:length(Edges)-1,'Uni',0));
ValPerGroupstd2 = cell2mat(arrayfun(@(X) nanstd(PercStableConservative(2,(PercStableConservative(1,:)>=Edges(X)&PercStableConservative(1,:)<=Edges(X+1)))),1:length(Edges)-1,'Uni',0));


h=barwitherr([ValPerGroupstd;ValPerGroupstd2],[Vals],[ValPerGroupMean;ValPerGroupMean2]);
ylabel('Fraction tracked neurons')
xlabel('\Delta days')
title('Tracking across days')
legend({'Liberal','Conservative'})
makepretty
offsetAxes


%% as a function of drift?
subplot(2,2,3)
scatter(PercStableLiberal(3,:),PercStableLiberal(2,:),20,[0 0 0],'filled')
hold on
scatter(PercStableConservative(3,:),PercStableConservative(2,:),20,[1 0 0],'filled')
ylabel('Fraction tracked neurons')
title('Tracking across drift')

xlabel('Drift (\muM)')
makepretty
offsetAxes

subplot(2,2,4)
hold on
Edges = min(PercStableLiberal(3,:)):5:max(PercStableLiberal(3,:));
Vals =  min(PercStableLiberal(3,:))+2.5:5:max(PercStableLiberal(3,:))-2.5;
ValPerGroupMean = cell2mat(arrayfun(@(X) nanmean(PercStableLiberal(2,(PercStableLiberal(3,:)>=Edges(X)&PercStableLiberal(3,:)<=Edges(X+1)))),1:length(Edges)-1,'Uni',0));
ValPerGroupstd = cell2mat(arrayfun(@(X) nanstd(PercStableLiberal(2,(PercStableLiberal(3,:)>=Edges(X)&PercStableLiberal(3,:)<=Edges(X+1)))),1:length(Edges)-1,'Uni',0));
ValPerGroupMean2 = cell2mat(arrayfun(@(X) nanmean(PercStableConservative(2,(PercStableConservative(3,:)>=Edges(X)&PercStableConservative(3,:)<=Edges(X+1)))),1:length(Edges)-1,'Uni',0));
ValPerGroupstd2 = cell2mat(arrayfun(@(X) nanstd(PercStableConservative(2,(PercStableConservative(3,:)>=Edges(X)&PercStableConservative(3,:)<=Edges(X+1)))),1:length(Edges)-1,'Uni',0));


h=barwitherr([ValPerGroupstd;ValPerGroupstd2],[Vals],[ValPerGroupMean;ValPerGroupMean2]);
ylabel('Fraction tracked neurons')
xlabel('Drift (\muM)')
title('Tracking across drift')
legend('Liberal','Conservative')
makepretty
offsetAxes



save(fullfile(UMparam.SaveDir,'PercentageStable.mat'),'PercStableLiberal','PercStableConservative')


