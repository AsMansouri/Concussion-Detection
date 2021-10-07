%% EEG Analysis
close all;
clear;
clc;

addpath(genpath('..\Functions'));
addpath(genpath('..\Functions\FastICA_25'));

%% Loading Participants Info
load('PrePostInfos.mat');
PrePostParticipants([8,9,11]) = [];

for i = 1:length(PrePostParticipants)
        
    for j = 1:length(PrePostParticipants(i).Sessions)
        PrePostParticipants(i).Sessions(j).nbackPartition = replace(PrePostParticipants(i).Sessions(j).nbackPartition , ...
            'J:\Concussion\Dataset\TBIC edf- baselinecorrected_Clean_Partitioned\', 'Partitions/'); 
        load(PrePostParticipants(i).Sessions(j).nbackPartition);
        PrePostParticipants(i).Sessions(j).FreqBand = FreqBand;
        clear FreqBand;
    end
    
end

% Calculating Po2Corr3 and Po2Corr2 distances
%%
parfor i = 1:length(PrePostParticipants)
    disp(num2str(PrePostParticipants(i).ID)); 
    for f = 1:length(PrePostParticipants(i).Sessions(1).FreqBand)
      
        for j = 1:length(PrePostParticipants(i).Sessions)
        
            tmpPower = sum(PrePostParticipants(i).Sessions(j).FreqBand(f).fftDataset.^2,3);
                          
            tmpEuclid = PrePostParticipants(i).Sessions(j).FreqBand(f).EuclidDistDatasets;
            tmpEuclid = (tmpEuclid - min(tmpEuclid(:))) / (max(tmpEuclid(:)) - min(tmpEuclid(:))); % normalizing Euclid distances!!
            
            tmpCos = PrePostParticipants(i).Sessions(j).FreqBand(f).CosDistDatasets;
            tmpCos = (tmpCos - min(tmpCos(:))) / (max(tmpCos(:)) - min(tmpCos(:))); % normalizing Cos distances!!
            
            tmpCorr = PrePostParticipants(i).Sessions(j).FreqBand(f).CorrDistDatasets;
            tmpCorr = (tmpCorr - min(tmpCorr(:))) / (max(tmpCorr(:)) - min(tmpCorr(:))); % normalizing Corr distances!!
            
%             tmpPrePo2Euclid3 = zeros(size(tmpEuclid,1),size(tmpPower,2),size(tmpPower,2));
%             tmpPrePo2Euclid2 = zeros(size(tmpEuclid,1),size(tmpPower,2),size(tmpPower,2));
%             tmpPrePo3Euclid2 = zeros(size(tmpEuclid,1),size(tmpPower,2),size(tmpPower,2));
%             
%             tmpPrePo2Cos3 = zeros(size(tmpCos,1),size(tmpPower,2),size(tmpPower,2));
%             tmpPrePo2Cos2 = zeros(size(tmpCos,1),size(tmpPower,2),size(tmpPower,2));
%             tmpPrePo3Cos2 = zeros(size(tmpCos,1),size(tmpPower,2),size(tmpPower,2));
%             
%             tmpPrePo2Corr3 = zeros(size(tmpCorr,1),size(tmpPower,2),size(tmpPower,2));
%             tmpPrePo2Corr2 = zeros(size(tmpCorr,1),size(tmpPower,2),size(tmpPower,2));
%             tmpPrePo3Corr2 = zeros(size(tmpCorr,1),size(tmpPower,2),size(tmpPower,2));
            
            tmpAll = nthroot(tmpEuclid.^3 + tmpCos.^3 + tmpCorr.^3,3);            
            tmpAllwPo = zeros(size(tmpEuclid,1),size(tmpPower,2),size(tmpPower,2));
            
            for seg = 1:size(tmpCorr,1)
                tmpnormPower = [];
                tmpnormPower(1,:) = tmpPower(seg,:)/max(tmpPower(seg,:));  % Normalizng POwers
                for p1 = 1:size(tmpPower,2)
                    for p2 = 1:size(tmpPower,2)
%                         tmpPrePo2Euclid3(seg,p1,p2) = sqrt((tmpEuclid(seg,p1,p2)).^3 + mean([tmpnormPower(p1) tmpnormPower(p2)]).^2);
%                         tmpPrePo2Euclid2(seg,p1,p2) = sqrt((tmpEuclid(seg,p1,p2)).^2 + mean([tmpnormPower(p1) tmpnormPower(p2)]).^2);  
%                         tmpPrePo3Euclid2(seg,p1,p2) = sqrt((tmpEuclid(seg,p1,p2)).^2 + mean([tmpnormPower(p1) tmpnormPower(p2)]).^3);  
%                         
%                         tmpPrePo2Cos3(seg,p1,p2) = sqrt((tmpCos(seg,p1,p2)).^3 + mean([tmpnormPower(p1) tmpnormPower(p2)]).^2);
%                         tmpPrePo2Cos2(seg,p1,p2) = sqrt((tmpCos(seg,p1,p2)).^2 + mean([tmpnormPower(p1) tmpnormPower(p2)]).^2);  
%                         tmpPrePo3Cos2(seg,p1,p2) = sqrt((tmpCos(seg,p1,p2)).^2 + mean([tmpnormPower(p1) tmpnormPower(p2)]).^3); 
%                         
%                         tmpPrePo2Corr3(seg,p1,p2) = sqrt((tmpCorr(seg,p1,p2)).^3 + mean([tmpnormPower(p1) tmpnormPower(p2)]).^2);
%                         tmpPrePo2Corr2(seg,p1,p2) = sqrt((tmpCorr(seg,p1,p2)).^2 + mean([tmpnormPower(p1) tmpnormPower(p2)]).^2);  
%                         tmpPrePo3Corr2(seg,p1,p2) = sqrt((tmpCorr(seg,p1,p2)).^2 + mean([tmpnormPower(p1) tmpnormPower(p2)]).^3); 
                        
                        tmpAllwPo(seg,p1,p2) = nthroot(( (tmpEuclid(seg,p1,p2)).^4 + (tmpCos(seg,p1,p2)).^4 + ...
                                            (tmpCorr(seg,p1,p2)).^4 + mean([tmpnormPower(p1) tmpnormPower(p2)]).^4), 4); 
                        
                    end
                end
            end
%             tmpPrePo2Euclid3 = tmpPrePo2Euclid3 / max(tmpPrePo2Euclid3(:));
%             tmpPrePo2Euclid2 = tmpPrePo2Euclid2 / max(tmpPrePo2Euclid2(:));
%             tmpPrePo3Euclid2 = tmpPrePo3Euclid2 / max(tmpPrePo3Euclid2(:));
%             
%             PrePostParticipants(i).Sessions(j).FreqBand(f).Po2Euclid3DistDatasets = tmpPrePo2Euclid3;
%             PrePostParticipants(i).Sessions(j).FreqBand(f).Po2Euclid2DistDatasets = tmpPrePo2Euclid2; 
%             PrePostParticipants(i).Sessions(j).FreqBand(f).Po3Euclid2DistDatasets = tmpPrePo3Euclid2; 
%             
%             tmpPrePo2Cos3 = tmpPrePo2Cos3 / max(tmpPrePo2Cos3(:));
%             tmpPrePo2Cos2 = tmpPrePo2Cos2 / max(tmpPrePo2Cos2(:));
%             tmpPrePo3Cos2 = tmpPrePo3Cos2 / max(tmpPrePo3Cos2(:));
% 
%             PrePostParticipants(i).Sessions(j).FreqBand(f).Po2Cos3DistDatasets = tmpPrePo2Cos3;
%             PrePostParticipants(i).Sessions(j).FreqBand(f).Po2Cos2DistDatasets = tmpPrePo2Cos2; 
%             PrePostParticipants(i).Sessions(j).FreqBand(f).Po3Cos2DistDatasets = tmpPrePo3Cos2; 
%             
%             tmpPrePo2Corr3 = tmpPrePo2Corr3 / max(tmpPrePo2Corr3(:));
%             tmpPrePo2Corr2 = tmpPrePo2Corr2 / max(tmpPrePo2Corr2(:));
%             tmpPrePo3Corr2 = tmpPrePo3Corr2 / max(tmpPrePo3Corr2(:));
% 
%             PrePostParticipants(i).Sessions(j).FreqBand(f).Po2Corr3DistDatasets = tmpPrePo2Corr3;
%             PrePostParticipants(i).Sessions(j).FreqBand(f).Po2Corr2DistDatasets = tmpPrePo2Corr2; 
%             PrePostParticipants(i).Sessions(j).FreqBand(f).Po3Corr2DistDatasets = tmpPrePo3Corr2; 
            
            PrePostParticipants(i).Sessions(j).FreqBand(f).AllFeatDistDatasets = tmpAll;
            PrePostParticipants(i).Sessions(j).FreqBand(f).AllFeatwPoDistDatasets = tmpAllwPo; 
        end
    end
end
        

%% Generating Datasets
for dsnum = 1:10
    mkdir(['Dataset',num2str(dsnum),'/Datasets_UnRanked/']);
    mkdir(['Dataset',num2str(dsnum),'/Datasets_TimeRanked/']);
    mkdir(['Dataset',num2str(dsnum),'/Datasets_FreqRanked/']);
    mkdir(['Dataset',num2str(dsnum),'/Datasets_TimeFreqRanked/']);
    
    disp(['Dataset',num2str(dsnum),':']);
    for i = 1:length(PrePostParticipants)
        disp(num2str(PrePostParticipants(i).ID)); 
        tic;
        for f = 1:length(PrePostParticipants(i).Sessions(1).FreqBand)

            %% UnRanked

            % - AllFeatures
            Dataset_Train=[];Labels_Train=[];Dataset_Test=[];Labels_Test=[];

            tmpPre = PrePostParticipants(i).Sessions(1).FreqBand(f).AllFeatDistDatasets;
            tmpPre(PrePostParticipants(i).Sessions(1).BadObsSegs,:,:) = [];
            tmpPost = PrePostParticipants(i).Sessions(2).FreqBand(f).AllFeatDistDatasets;
            tmpPost(PrePostParticipants(i).Sessions(2).BadObsSegs,:,:) = [];
            NumSeg = min(size(tmpPre,1),size(tmpPost,1));

            tmpPermInd = randperm(size(tmpPre,1));
            Pre = tmpPre(tmpPermInd(1:NumSeg),:,:);
            tmpPermInd = randperm(size(tmpPost,1));
            Post = tmpPost(tmpPermInd(1:NumSeg),:,:);

            Dataset_Test = [Pre,Post];

            if PrePostParticipants(i).Label == 1
                Labels_Test = ones(NumSeg,1);
            else
                Labels_Test = zeros(NumSeg,1);
            end

            for j = 1:length(PrePostParticipants)
                if j == i 
                    continue;
                end

                tmpPre = PrePostParticipants(j).Sessions(1).FreqBand(f).AllFeatDistDatasets;
                tmpPre(PrePostParticipants(j).Sessions(1).BadObsSegs,:,:) = [];
                tmpPost = PrePostParticipants(j).Sessions(2).FreqBand(f).AllFeatDistDatasets;
                tmpPost(PrePostParticipants(j).Sessions(2).BadObsSegs,:,:) = [];
                NumSeg = min(size(tmpPre,1),size(tmpPost,1));

                tmpPermInd = randperm(size(tmpPre,1));
                Pre = tmpPre(tmpPermInd(1:NumSeg),:,:);
                tmpPermInd = randperm(size(tmpPost,1));
                Post = tmpPost(tmpPermInd(1:NumSeg),:,:);

                tmpTrain = [Pre,Post];

                Dataset_Train = [Dataset_Train; tmpTrain];

                if PrePostParticipants(j).Label == 1
                    tmpLabels = ones(NumSeg,1);
                else
                    tmpLabels = zeros(NumSeg,1);
                end

                Labels_Train = [Labels_Train;tmpLabels];

            end

            try
                save(['Dataset',num2str(dsnum),'/Datasets_UnRanked/Datasets_CleanSegments_PrePost_Freq_AllFeat_',PrePostParticipants(i).Sessions(1).FreqBand(f).Name,'_ID_',num2str(PrePostParticipants(i).ID),'.mat'],...
                      'Dataset_Train','Dataset_Test','Labels_Train','Labels_Test');
            catch
                warning('Saved as version v7.3......');
                save(['Dataset',num2str(dsnum),'/Datasets_UnRanked/Datasets_CleanSegments_PrePost_Freq_AllFeat_',PrePostParticipants(i).Sessions(1).FreqBand(f).Name,'_ID_',num2str(PrePostParticipants(i).ID),'V73.mat'],...
                      'Dataset_Train','Dataset_Test','Labels_Train','Labels_Test','-v7.3');
            end


            % - AllFeatures with Power
            Dataset_Train=[];Labels_Train=[];Dataset_Test=[];Labels_Test=[];

            tmpPre = PrePostParticipants(i).Sessions(1).FreqBand(f).AllFeatwPoDistDatasets;
            tmpPre(PrePostParticipants(i).Sessions(1).BadObsSegs,:,:) = [];
            tmpPost = PrePostParticipants(i).Sessions(2).FreqBand(f).AllFeatwPoDistDatasets;
            tmpPost(PrePostParticipants(i).Sessions(2).BadObsSegs,:,:) = [];
            NumSeg = min(size(tmpPre,1),size(tmpPost,1));

            tmpPermInd = randperm(size(tmpPre,1));
            Pre = tmpPre(tmpPermInd(1:NumSeg),:,:);
            tmpPermInd = randperm(size(tmpPost,1));
            Post = tmpPost(tmpPermInd(1:NumSeg),:,:);

            Dataset_Test = [Pre,Post];

            if PrePostParticipants(i).Label == 1
                Labels_Test = ones(NumSeg,1);
            else
                Labels_Test = zeros(NumSeg,1);
            end

            for j = 1:length(PrePostParticipants)
                if j == i 
                    continue;
                end

                tmpPre = PrePostParticipants(j).Sessions(1).FreqBand(f).AllFeatwPoDistDatasets;
                tmpPre(PrePostParticipants(j).Sessions(1).BadObsSegs,:,:) = [];
                tmpPost = PrePostParticipants(j).Sessions(2).FreqBand(f).AllFeatwPoDistDatasets;
                tmpPost(PrePostParticipants(j).Sessions(2).BadObsSegs,:,:) = [];
                NumSeg = min(size(tmpPre,1),size(tmpPost,1));

                tmpPermInd = randperm(size(tmpPre,1));
                Pre = tmpPre(tmpPermInd(1:NumSeg),:,:);
                tmpPermInd = randperm(size(tmpPost,1));
                Post = tmpPost(tmpPermInd(1:NumSeg),:,:);

                tmpTrain = [Pre,Post];

                Dataset_Train = [Dataset_Train; tmpTrain];

                if PrePostParticipants(j).Label == 1
                    tmpLabels = ones(NumSeg,1);
                else
                    tmpLabels = zeros(NumSeg,1);
                end

                Labels_Train = [Labels_Train;tmpLabels];

            end

            try
                save(['Dataset',num2str(dsnum),'/Datasets_UnRanked/Datasets_CleanSegments_PrePost_Freq_AllwPo_',PrePostParticipants(i).Sessions(1).FreqBand(f).Name,'_ID_',num2str(PrePostParticipants(i).ID),'.mat'],...
                      'Dataset_Train','Dataset_Test','Labels_Train','Labels_Test');
            catch
                warning('Saved as version v7.3......');
                save(['Dataset',num2str(dsnum),'/Datasets_UnRanked/Datasets_CleanSegments_PrePost_Freq_AllwPo_',PrePostParticipants(i).Sessions(1).FreqBand(f).Name,'_ID_',num2str(PrePostParticipants(i).ID),'V73.mat'],...
                      'Dataset_Train','Dataset_Test','Labels_Train','Labels_Test','-v7.3');
            end

            %% Time Ranked

            % - AllFeatures
            Dataset_Train=[];Labels_Train=[];Dataset_Test=[];Labels_Test=[];

            tmpPre = PrePostParticipants(i).Sessions(1).FreqBand(f).AllFeatDistDatasets;
            tmpPost = PrePostParticipants(i).Sessions(2).FreqBand(f).AllFeatDistDatasets;
            NumSeg = 45;

            tmpPermInd = randperm(NumSeg);
            Pre = tmpPre((PrePostParticipants(i).Sessions(1).TimeRank1(tmpPermInd)),:,:);
            tmpPermInd = randperm(NumSeg);
            Post = tmpPost((PrePostParticipants(i).Sessions(2).TimeRank1(tmpPermInd)),:,:);

            Dataset_Test = [Pre,Post];

            if PrePostParticipants(i).Label == 1
                Labels_Test = ones(NumSeg,1);
            else
                Labels_Test = zeros(NumSeg,1);
            end

            for j = 1:length(PrePostParticipants)
                if j == i 
                    continue;
                end

                tmpPre = PrePostParticipants(j).Sessions(1).FreqBand(f).AllFeatDistDatasets;
                tmpPost = PrePostParticipants(j).Sessions(2).FreqBand(f).AllFeatDistDatasets;
                NumSeg = 45;

                tmpPermInd = randperm(NumSeg);
                Pre = tmpPre((PrePostParticipants(j).Sessions(1).TimeRank1(tmpPermInd)),:,:);
                tmpPermInd = randperm(NumSeg);
                Post = tmpPost((PrePostParticipants(j).Sessions(2).TimeRank1(tmpPermInd)),:,:);

                tmpTrain = [Pre,Post];

                Dataset_Train = [Dataset_Train; tmpTrain];

                if PrePostParticipants(j).Label == 1
                    tmpLabels = ones(NumSeg,1);
                else
                    tmpLabels = zeros(NumSeg,1);
                end

                Labels_Train = [Labels_Train;tmpLabels];

            end

            try
                save(['Dataset',num2str(dsnum),'/Datasets_TimeRanked/Datasets_CleanSegments_PrePost_Freq_AllFeat_',PrePostParticipants(i).Sessions(1).FreqBand(f).Name,'_ID_',num2str(PrePostParticipants(i).ID),'.mat'],...
                      'Dataset_Train','Dataset_Test','Labels_Train','Labels_Test');
            catch
                warning('Saved as version v7.3......');
                save(['Dataset',num2str(dsnum),'/Datasets_TimeRanked/Datasets_CleanSegments_PrePost_Freq_AllFeat_',PrePostParticipants(i).Sessions(1).FreqBand(f).Name,'_ID_',num2str(PrePostParticipants(i).ID),'V73.mat'],...
                      'Dataset_Train','Dataset_Test','Labels_Train','Labels_Test','-v7.3');
            end


            % - AllFeatures with Power
            Dataset_Train=[];Labels_Train=[];Dataset_Test=[];Labels_Test=[];

            tmpPre = PrePostParticipants(i).Sessions(1).FreqBand(f).AllFeatwPoDistDatasets;
            tmpPost = PrePostParticipants(i).Sessions(2).FreqBand(f).AllFeatwPoDistDatasets;
            NumSeg = 45;

            tmpPermInd = randperm(NumSeg);
            Pre = tmpPre((PrePostParticipants(i).Sessions(1).TimeRank1(tmpPermInd)),:,:);
            tmpPermInd = randperm(NumSeg);
            Post = tmpPost((PrePostParticipants(i).Sessions(2).TimeRank1(tmpPermInd)),:,:);

            Dataset_Test = [Pre,Post];

            if PrePostParticipants(i).Label == 1
                Labels_Test = ones(NumSeg,1);
            else
                Labels_Test = zeros(NumSeg,1);
            end

            for j = 1:length(PrePostParticipants)
                if j == i 
                    continue;
                end

                tmpPre = PrePostParticipants(j).Sessions(1).FreqBand(f).AllFeatwPoDistDatasets;
                tmpPost = PrePostParticipants(j).Sessions(2).FreqBand(f).AllFeatwPoDistDatasets;
                NumSeg = 45;

                tmpPermInd = randperm(NumSeg);
                Pre = tmpPre((PrePostParticipants(j).Sessions(1).TimeRank1(tmpPermInd)),:,:);
                tmpPermInd = randperm(NumSeg);
                Post = tmpPost((PrePostParticipants(j).Sessions(2).TimeRank1(tmpPermInd)),:,:);

                tmpTrain = [Pre,Post];

                Dataset_Train = [Dataset_Train; tmpTrain];

                if PrePostParticipants(j).Label == 1
                    tmpLabels = ones(NumSeg,1);
                else
                    tmpLabels = zeros(NumSeg,1);
                end

                Labels_Train = [Labels_Train;tmpLabels];

            end

            try
                save(['Dataset',num2str(dsnum),'/Datasets_TimeRanked/Datasets_CleanSegments_PrePost_Freq_AllwPo_',PrePostParticipants(i).Sessions(1).FreqBand(f).Name,'_ID_',num2str(PrePostParticipants(i).ID),'.mat'],...
                      'Dataset_Train','Dataset_Test','Labels_Train','Labels_Test');
            catch
                warning('Saved as version v7.3......');
                save(['Dataset',num2str(dsnum),'/Datasets_TimeRanked/Datasets_CleanSegments_PrePost_Freq_AllwPo_',PrePostParticipants(i).Sessions(1).FreqBand(f).Name,'_ID_',num2str(PrePostParticipants(i).ID),'V73.mat'],...
                      'Dataset_Train','Dataset_Test','Labels_Train','Labels_Test','-v7.3');
            end

            %% Freq Ranked

            % - AllFeatures
            Dataset_Train=[];Labels_Train=[];Dataset_Test=[];Labels_Test=[];

            tmpPre = PrePostParticipants(i).Sessions(1).FreqBand(f).AllFeatDistDatasets;
            tmpPost = PrePostParticipants(i).Sessions(2).FreqBand(f).AllFeatDistDatasets;
            NumSeg = 45;

            tmpPermInd = randperm(NumSeg);
            Pre = tmpPre((PrePostParticipants(i).Sessions(1).FreqRank1(tmpPermInd)),:,:);
            tmpPermInd = randperm(NumSeg);
            Post = tmpPost((PrePostParticipants(i).Sessions(2).FreqRank1(tmpPermInd)),:,:);

            Dataset_Test = [Pre,Post];

            if PrePostParticipants(i).Label == 1
                Labels_Test = ones(NumSeg,1);
            else
                Labels_Test = zeros(NumSeg,1);
            end

            for j = 1:length(PrePostParticipants)
                if j == i 
                    continue;
                end

                tmpPre = PrePostParticipants(j).Sessions(1).FreqBand(f).AllFeatDistDatasets;
                tmpPost = PrePostParticipants(j).Sessions(2).FreqBand(f).AllFeatDistDatasets;
                NumSeg = 45;

                tmpPermInd = randperm(NumSeg);
                Pre = tmpPre((PrePostParticipants(j).Sessions(1).FreqRank1(tmpPermInd)),:,:);
                tmpPermInd = randperm(NumSeg);
                Post = tmpPost((PrePostParticipants(j).Sessions(2).FreqRank1(tmpPermInd)),:,:);

                tmpTrain = [Pre,Post];

                Dataset_Train = [Dataset_Train; tmpTrain];

                if PrePostParticipants(j).Label == 1
                    tmpLabels = ones(NumSeg,1);
                else
                    tmpLabels = zeros(NumSeg,1);
                end

                Labels_Train = [Labels_Train;tmpLabels];

            end

            try
                save(['Dataset',num2str(dsnum),'/Datasets_FreqRanked/Datasets_CleanSegments_PrePost_Freq_AllFeat_',PrePostParticipants(i).Sessions(1).FreqBand(f).Name,'_ID_',num2str(PrePostParticipants(i).ID),'.mat'],...
                      'Dataset_Train','Dataset_Test','Labels_Train','Labels_Test');
            catch
                warning('Saved as version v7.3......');
                save(['Dataset',num2str(dsnum),'/Datasets_FreqRanked/Datasets_CleanSegments_PrePost_Freq_AllFeat_',PrePostParticipants(i).Sessions(1).FreqBand(f).Name,'_ID_',num2str(PrePostParticipants(i).ID),'V73.mat'],...
                      'Dataset_Train','Dataset_Test','Labels_Train','Labels_Test','-v7.3');
            end


            % - AllFeatures with Power
            Dataset_Train=[];Labels_Train=[];Dataset_Test=[];Labels_Test=[];

            tmpPre = PrePostParticipants(i).Sessions(1).FreqBand(f).AllFeatwPoDistDatasets;
            tmpPost = PrePostParticipants(i).Sessions(2).FreqBand(f).AllFeatwPoDistDatasets;
            NumSeg = 45;

            tmpPermInd = randperm(NumSeg);
            Pre = tmpPre((PrePostParticipants(i).Sessions(1).FreqRank1(tmpPermInd)),:,:);
            tmpPermInd = randperm(NumSeg);
            Post = tmpPost((PrePostParticipants(i).Sessions(2).FreqRank1(tmpPermInd)),:,:);

            Dataset_Test = [Pre,Post];

            if PrePostParticipants(i).Label == 1
                Labels_Test = ones(NumSeg,1);
            else
                Labels_Test = zeros(NumSeg,1);
            end

            for j = 1:length(PrePostParticipants)
                if j == i 
                    continue;
                end

                tmpPre = PrePostParticipants(j).Sessions(1).FreqBand(f).AllFeatwPoDistDatasets;
                tmpPost = PrePostParticipants(j).Sessions(2).FreqBand(f).AllFeatwPoDistDatasets;
                NumSeg = 45;

                tmpPermInd = randperm(NumSeg);
                Pre = tmpPre((PrePostParticipants(j).Sessions(1).FreqRank1(tmpPermInd)),:,:);
                tmpPermInd = randperm(NumSeg);
                Post = tmpPost((PrePostParticipants(j).Sessions(2).FreqRank1(tmpPermInd)),:,:);

                tmpTrain = [Pre,Post];

                Dataset_Train = [Dataset_Train; tmpTrain];

                if PrePostParticipants(j).Label == 1
                    tmpLabels = ones(NumSeg,1);
                else
                    tmpLabels = zeros(NumSeg,1);
                end

                Labels_Train = [Labels_Train;tmpLabels];

            end

            try
                save(['Dataset',num2str(dsnum),'/Datasets_FreqRanked/Datasets_CleanSegments_PrePost_Freq_AllwPo_',PrePostParticipants(i).Sessions(1).FreqBand(f).Name,'_ID_',num2str(PrePostParticipants(i).ID),'.mat'],...
                      'Dataset_Train','Dataset_Test','Labels_Train','Labels_Test');
            catch
                warning('Saved as version v7.3......');
                save(['Dataset',num2str(dsnum),'/Datasets_FreqRanked/Datasets_CleanSegments_PrePost_Freq_AllwPo_',PrePostParticipants(i).Sessions(1).FreqBand(f).Name,'_ID_',num2str(PrePostParticipants(i).ID),'V73.mat'],...
                      'Dataset_Train','Dataset_Test','Labels_Train','Labels_Test','-v7.3');
            end

            %% TimeFreq Ranked

            % - AllFeatures
            Dataset_Train=[];Labels_Train=[];Dataset_Test=[];Labels_Test=[];

            tmpPre = PrePostParticipants(i).Sessions(1).FreqBand(f).AllFeatDistDatasets;
            tmpPost = PrePostParticipants(i).Sessions(2).FreqBand(f).AllFeatDistDatasets;
            NumSeg = 45;

            tmpPermInd = randperm(NumSeg);
            Pre = tmpPre((PrePostParticipants(i).Sessions(1).RankFreqTime1(tmpPermInd)),:,:);
            tmpPermInd = randperm(NumSeg);
            Post = tmpPost((PrePostParticipants(i).Sessions(2).RankFreqTime1(tmpPermInd)),:,:);

            Dataset_Test = [Pre,Post];

            if PrePostParticipants(i).Label == 1
                Labels_Test = ones(NumSeg,1);
            else
                Labels_Test = zeros(NumSeg,1);
            end

            for j = 1:length(PrePostParticipants)
                if j == i 
                    continue;
                end

                tmpPre = PrePostParticipants(j).Sessions(1).FreqBand(f).AllFeatDistDatasets;
                tmpPost = PrePostParticipants(j).Sessions(2).FreqBand(f).AllFeatDistDatasets;
                NumSeg = 45;

                tmpPermInd = randperm(NumSeg);
                Pre = tmpPre((PrePostParticipants(j).Sessions(1).RankFreqTime1(tmpPermInd)),:,:);
                tmpPermInd = randperm(NumSeg);
                Post = tmpPost((PrePostParticipants(j).Sessions(2).RankFreqTime1(tmpPermInd)),:,:);

                tmpTrain = [Pre,Post];

                Dataset_Train = [Dataset_Train; tmpTrain];

                if PrePostParticipants(j).Label == 1
                    tmpLabels = ones(NumSeg,1);
                else
                    tmpLabels = zeros(NumSeg,1);
                end

                Labels_Train = [Labels_Train;tmpLabels];

            end

            try
                save(['Dataset',num2str(dsnum),'/Datasets_TimeFreqRanked/Datasets_CleanSegments_PrePost_Freq_AllFeat_',PrePostParticipants(i).Sessions(1).FreqBand(f).Name,'_ID_',num2str(PrePostParticipants(i).ID),'.mat'],...
                      'Dataset_Train','Dataset_Test','Labels_Train','Labels_Test');
            catch
                warning('Saved as version v7.3......');
                save(['Dataset',num2str(dsnum),'/Datasets_TimeFreqRanked/Datasets_CleanSegments_PrePost_Freq_AllFeat_',PrePostParticipants(i).Sessions(1).FreqBand(f).Name,'_ID_',num2str(PrePostParticipants(i).ID),'V73.mat'],...
                      'Dataset_Train','Dataset_Test','Labels_Train','Labels_Test','-v7.3');
            end


            % - AllFeatures with Power
            Dataset_Train=[];Labels_Train=[];Dataset_Test=[];Labels_Test=[];

            tmpPre = PrePostParticipants(i).Sessions(1).FreqBand(f).AllFeatwPoDistDatasets;
            tmpPost = PrePostParticipants(i).Sessions(2).FreqBand(f).AllFeatwPoDistDatasets;
            NumSeg = 45;

            tmpPermInd = randperm(NumSeg);
            Pre = tmpPre((PrePostParticipants(i).Sessions(1).RankFreqTime1(tmpPermInd)),:,:);
            tmpPermInd = randperm(NumSeg);
            Post = tmpPost((PrePostParticipants(i).Sessions(2).RankFreqTime1(tmpPermInd)),:,:);

            Dataset_Test = [Pre,Post];

            if PrePostParticipants(i).Label == 1
                Labels_Test = ones(NumSeg,1);
            else
                Labels_Test = zeros(NumSeg,1);
            end

            for j = 1:length(PrePostParticipants)
                if j == i 
                    continue;
                end

                tmpPre = PrePostParticipants(j).Sessions(1).FreqBand(f).AllFeatwPoDistDatasets;
                tmpPost = PrePostParticipants(j).Sessions(2).FreqBand(f).AllFeatwPoDistDatasets;
                NumSeg = 45;

                tmpPermInd = randperm(NumSeg);
                Pre = tmpPre((PrePostParticipants(j).Sessions(1).RankFreqTime1(tmpPermInd)),:,:);
                tmpPermInd = randperm(NumSeg);
                Post = tmpPost((PrePostParticipants(j).Sessions(2).RankFreqTime1(tmpPermInd)),:,:);

                tmpTrain = [Pre,Post];

                Dataset_Train = [Dataset_Train; tmpTrain];

                if PrePostParticipants(j).Label == 1
                    tmpLabels = ones(NumSeg,1);
                else
                    tmpLabels = zeros(NumSeg,1);
                end

                Labels_Train = [Labels_Train;tmpLabels];

            end

            try
                save(['Dataset',num2str(dsnum),'/Datasets_TimeFreqRanked/Datasets_CleanSegments_PrePost_Freq_AllwPo_',PrePostParticipants(i).Sessions(1).FreqBand(f).Name,'_ID_',num2str(PrePostParticipants(i).ID),'.mat'],...
                      'Dataset_Train','Dataset_Test','Labels_Train','Labels_Test');
            catch
                warning('Saved as version v7.3......');
                save(['Dataset',num2str(dsnum),'/Datasets_TimeFreqRanked/Datasets_CleanSegments_PrePost_Freq_AllwPo_',PrePostParticipants(i).Sessions(1).FreqBand(f).Name,'_ID_',num2str(PrePostParticipants(i).ID),'V73.mat'],...
                      'Dataset_Train','Dataset_Test','Labels_Train','Labels_Test','-v7.3');
            end


        end
        toc
    end
end        

%%
disp('-----Done!!-----');



%% Functions
function [DistData] = DistGenerator(Data, DistMetricType)

    tmp = size(Data);
    DistData = zeros(tmp(1),tmp(2),tmp(2));

    for seg = 1:size(Data,1)
        tmp11(:,:) = Data(seg,:,:);
        D = pdist(tmp11, DistMetricType);
        DistData(seg,:,:) = squareform(D);
    end
end

function [AMIs] = AMIElec(Data,NLevel,MAXAMIK)
%[AMIs] = AMI(Data,NLevel,MAXAMIK)
%   Calculate the AMI by uniformly quantization. The Data should be a row
%   vector or a row matrix.
    
    Mins = min(Data(:));
    Maxs = max(Data(:));
    
    NormData = ((2.*Data)-(Maxs+Mins))./(Maxs-Mins);
    Data = NormData;
    Levels = [-1:((1)-(-1))/NLevel:1];
    
    NRow = size(Data,1);
    TotalSamp = size(Data,2);
    AMIs = zeros(NRow, MAXAMIK);
    JointProb = zeros(NLevel);
    Prob = zeros(1,NLevel);

    for m = 1:NLevel
       tmp = Data >= Levels(m) & Data < Levels(m+1);
       Prob(m) = sum(tmp(:))/length(Data(:));
    end

    Loc = cell(1,NLevel);
    for i = 1:NRow
        for m = 1:NLevel
            Loc{m} = find(Data(i,:) >= Levels(m) & Data(i,:) < Levels(m+1));
        end
        for k = 1:MAXAMIK
            for p = 1:NLevel
                for q = 1:NLevel
                    [Lia,~] = ismember(Loc{p},Loc{q}-k);
                    JointProb(p,q) = sum(Lia)/(TotalSamp-k);
                    if Prob(p)~= 0 && Prob(q)~=0 && JointProb(p,q)~=0
                        AMIs(i,k) = AMIs(i,k) + JointProb(p,q)*log2(JointProb(p,q)/(Prob(p)*Prob(q)));
                    end
                end
            end
        end
    end
end

function [Dataset_Train,Labels_Train,Dataset_Test,Labels_Test] = DatasetGeneratorUnbal(Dataset, Labels)
    %Labels are 0 or 1
    Labels_Train = Labels;
    Dataset_Train = Dataset;
    Dataset_Test = []; Labels_Test = [];
    
    tmp1 = find(Labels_Train==1);
    tmp11 = randsample(tmp1,floor(length(tmp1)*.2));

    tmp0 = find(Labels_Train==0);
    tmp00 = randsample(tmp0,floor(length(tmp0)*.2));

    for t = 1:length(tmp11)
        Dataset_Test(end+1,:,:) = Dataset_Train(tmp11(t),:,:);
        Labels_Test(end+1,1) = Labels_Train(tmp11(t));
    end
    for t = 1:length(tmp00)
        Dataset_Test(end+1,:,:) = Dataset_Train(tmp00(t),:,:);
        Labels_Test(end+1,1) = Labels_Train(tmp00(t));
    end
    tmp01 = sort([tmp00;tmp11],'descend');
    for t = 1:length(tmp01)
        Dataset_Train(tmp01(t),:,:) = [];
        Labels_Train(tmp01(t)) = [];
    end

end

function [Dataset_Train,Labels_Train,Dataset_Test,Labels_Test] = DatasetGeneratorBal(Dataset, Labels)
    %Labels are 0 or 1
    Labels_Train = Labels;
    Dataset_Train = Dataset;
    Dataset_Test = []; Labels_Test = [];

    tmp1 = find(Labels_Train==1);
    tmp11 = randsample(tmp1,floor(length(tmp1)*.2));

    tmp0 = find(Labels_Train==0);
    tmp00 = randsample(tmp0,floor(length(tmp1)*.2));

    for t = 1:length(tmp11)
        Dataset_Test(end+1,:,:) = Dataset_Train(tmp11(t),:,:);
        Labels_Test(end+1,1) = Labels_Train(tmp11(t));
    end
    for t = 1:length(tmp00)
        Dataset_Test(end+1,:,:) = Dataset_Train(tmp00(t),:,:);
        Labels_Test(end+1,1) = Labels_Train(tmp00(t));
    end

    tmp01 = sort([tmp00;tmp11],'descend');
    for t = 1:length(tmp01)
        Dataset_Train(tmp01(t),:,:) = [];
        Labels_Train(tmp01(t)) = [];
    end

    tmp1 = find(Labels_Train==1);
    tmp00 = randsample(find(Labels_Train==0),length(find(Labels_Train==0))-length(tmp1));
    Dataset_Train(tmp00,:,:) = [];
    Labels_Train(tmp00) = [];
end

function [SegmentedSampelsInfo, SegmentedSampelsInfoLists] = SegmentedFilesInfo(Samplelisting)

SegmentedSampelsInfo = {};
for i = 3:length(Samplelisting)
    SegmentedSampelsInfo{i-2}.Name = Samplelisting(i).name;
    SegmentedSampelsInfo{i-2}.Address = [Samplelisting(i).folder '\' Samplelisting(i).name];
    if SegmentedSampelsInfo{i-2}.Name(1)== 'R'
        SegmentedSampelsInfo{i-2}.Name(1:2) = '';
    end
    switch SegmentedSampelsInfo{i-2}.Name(5)
        case 'F'
            SegmentedSampelsInfo{i-2}.Sport = 'football';
        case 'L'
            SegmentedSampelsInfo{i-2}.Sport = 'Lacrosse';
        case 'S'
            SegmentedSampelsInfo{i-2}.Sport = 'soccer';
        case 'R'
            SegmentedSampelsInfo{i-2}.Sport = 'rugby';
        otherwise
            SegmentedSampelsInfo{i-2}.Sport = 'ND';
    end
    
    if SegmentedSampelsInfo{i-2}.Name(1:3) == 'AFP'
        SegmentedSampelsInfo{i-2}.IDType = 'AFP';
        SegmentedSampelsInfo{i-2}.ID = SegmentedSampelsInfo{i-2}.Name(5:9);
        switch SegmentedSampelsInfo{i-2}.Name(10)
            case 'f'
                SegmentedSampelsInfo{i-2}.Gender = 'female';
            case 'm'
                SegmentedSampelsInfo{i-2}.Gender = 'male';
        end
        tmpy = split(SegmentedSampelsInfo{i-2}.Name,'_');
        if length(tmpy) == 4
            SegmentedSampelsInfo{i-2}.Year = 1;
        else
            SegmentedSampelsInfo{i-2}.Year = str2double(tmpy{4});
        end
        SegmentedSampelsInfo{i-2}.Session = str2double(cell2mat(strtrim(extractBetween(SegmentedSampelsInfo{i-2}.Name,'ses','_'))));
        if SegmentedSampelsInfo{i-2}.Name(16) == '_'
            if SegmentedSampelsInfo{i-2}.Name(17) == 'o' || SegmentedSampelsInfo{i-2}.Name(17) == 'O'
                SegmentedSampelsInfo{i-2}.Task = 'oddball';
            elseif SegmentedSampelsInfo{i-2}.Name(17) == 'n' || SegmentedSampelsInfo{i-2}.Name(17) == 'N'
                SegmentedSampelsInfo{i-2}.Task = 'nback';
            else
                SegmentedSampelsInfo{i-2}.Task = 'ND';
            end
        else
            if contains(lower(SegmentedSampelsInfo{i-2}.Name),'bac')
                SegmentedSampelsInfo{i-2}.Task = 'nback';
            elseif contains(lower(SegmentedSampelsInfo{i-2}.Name),'odd')
                SegmentedSampelsInfo{i-2}.Task = 'oddball';
            else
                SegmentedSampelsInfo{i-2}.Task = 'both';
            end
        end
    else
        SegmentedSampelsInfo{i-2}.IDType = 'TBIC';
        SegmentedSampelsInfo{i-2}.ID = SegmentedSampelsInfo{i-2}.Name(7:11);
        switch SegmentedSampelsInfo{i-2}.Name(12)
            case 'f'
                SegmentedSampelsInfo{i-2}.Gender = 'female';
            case 'm'
                SegmentedSampelsInfo{i-2}.Gender = 'male';
        end
        tmpy = split(SegmentedSampelsInfo{i-2}.Name,'_');
        if length(tmpy) == 4
            SegmentedSampelsInfo{i-2}.Year = 1;
        else
            SegmentedSampelsInfo{i-2}.Year = str2double(tmpy{4});
        end
        SegmentedSampelsInfo{i-2}.Session = str2double(cell2mat(strtrim(extractBetween(SegmentedSampelsInfo{i-2}.Name,'ses','_'))));
        if SegmentedSampelsInfo{i-2}.Name(20) == '_'
            if SegmentedSampelsInfo{i-2}.Name(21) == 'o' || SegmentedSampelsInfo{i-2}.Name(21) == 'O'
                SegmentedSampelsInfo{i-2}.Task = 'oddball';
            elseif SegmentedSampelsInfo{i-2}.Name(21) == 'n' || SegmentedSampelsInfo{i-2}.Name(21) == 'N'
                SegmentedSampelsInfo{i-2}.Task = 'nback';
            else
                SegmentedSampelsInfo{i-2}.Task = 'ND';
            end
        else
            if contains(lower(SegmentedSampelsInfo{i-2}.Name),'bac')
                SegmentedSampelsInfo{i-2}.Task = 'nback';
            elseif contains(lower(SegmentedSampelsInfo{i-2}.Name),'odd')
                SegmentedSampelsInfo{i-2}.Task = 'oddball';
            else
                SegmentedSampelsInfo{i-2}.Task = 'both';
            end
        end
    end
end

SegmentedSampelsInfoLists = {};
for i = 3:length(Samplelisting)
    SegmentedSampelsInfoLists.Name{i-2} = Samplelisting(i).name;
    SegmentedSampelsInfoLists.Address{i-2} = [Samplelisting(i).folder '\' Samplelisting(i).name];
    if SegmentedSampelsInfoLists.Name{i-2}(1)== 'R'
        SegmentedSampelsInfoLists.Name{i-2}(1:2) = '';
    end
    switch SegmentedSampelsInfoLists.Name{i-2}(5)
        case 'F'
            SegmentedSampelsInfoLists.Sport{i-2} = 'football';
        case 'L'
            SegmentedSampelsInfoLists.Sport{i-2} = 'Lacrosse';
        case 'S'
            SegmentedSampelsInfoLists.Sport{i-2} = 'soccer';
        case 'R'
            SegmentedSampelsInfoLists.Sport{i-2} = 'rugby';
        otherwise
            SegmentedSampelsInfoLists.Sport{i-2} = 'ND';
    end
    
    if SegmentedSampelsInfoLists.Name{i-2}(1:3) == 'AFP'
        SegmentedSampelsInfoLists.IDType{i-2} = 'AFP';
        SegmentedSampelsInfoLists.ID{i-2} = SegmentedSampelsInfoLists.Name{i-2}(5:9);
        switch SegmentedSampelsInfoLists.Name{i-2}(10)
            case 'f'
                SegmentedSampelsInfoLists.Gender{i-2} = 'female';
            case 'm'
                SegmentedSampelsInfoLists.Gender{i-2} = 'male';
        end
        tmpy = split(SegmentedSampelsInfoLists.Name{i-2},'_');
        if length(tmpy) == 4
            SegmentedSampelsInfoLists.Year{i-2} = 1;
        else
            SegmentedSampelsInfoLists.Year{i-2} = str2double(tmpy{4});
        end
        SegmentedSampelsInfoLists.Session{i-2} = str2double(cell2mat(strtrim(extractBetween(SegmentedSampelsInfoLists.Name{i-2},'ses','_'))));
        if SegmentedSampelsInfoLists.Name{i-2}(16) == '_'
            if SegmentedSampelsInfoLists.Name{i-2}(17) == 'o' || SegmentedSampelsInfoLists.Name{i-2}(17) == 'O'
                SegmentedSampelsInfoLists.Task{i-2} = 'oddball';
            elseif SegmentedSampelsInfoLists.Name{i-2}(17) == 'n' || SegmentedSampelsInfoLists.Name{i-2}(17) == 'N'
                SegmentedSampelsInfoLists.Task{i-2} = 'nback';
            else
                SegmentedSampelsInfoLists.Task{i-2} = 'ND';
            end
        else
            if contains(lower(SegmentedSampelsInfo{i-2}.Name),'bac')
                SegmentedSampelsInfoLists.Task{i-2} = 'nback';
            elseif contains(lower(SegmentedSampelsInfo{i-2}.Name),'odd')
                SegmentedSampelsInfoLists.Task{i-2} = 'oddball';
            else
                SegmentedSampelsInfoLists.Task{i-2}= 'both';
            end
        end
    else
        SegmentedSampelsInfoLists.IDType{i-2} = 'TBIC';
        SegmentedSampelsInfoLists.ID{i-2} = SegmentedSampelsInfoLists.Name{i-2}(7:11);
        switch SegmentedSampelsInfoLists.Name{i-2}(12)
            case 'f'
                SegmentedSampelsInfoLists.Gender{i-2} = 'female';
            case 'm'
                SegmentedSampelsInfoLists.Gender{i-2} = 'male';
        end
                tmpy = split(SegmentedSampelsInfoLists.Name{i-2},'_');
        if length(tmpy) == 4
            SegmentedSampelsInfoLists.Year{i-2} = 1;
        else
            SegmentedSampelsInfoLists.Year{i-2} = str2double(tmpy{4});
        end
        SegmentedSampelsInfoLists.Session{i-2} = str2double(cell2mat(strtrim(extractBetween(SegmentedSampelsInfoLists.Name{i-2},'ses','_'))));
        if SegmentedSampelsInfoLists.Name{i-2}(20) == '_'
            if SegmentedSampelsInfoLists.Name{i-2}(21) == 'o' || SegmentedSampelsInfoLists.Name{i-2}(21) == 'O'
                SegmentedSampelsInfoLists.Task{i-2} = 'oddball';
            elseif SegmentedSampelsInfoLists.Name{i-2}(21) == 'n' || SegmentedSampelsInfoLists.Name{i-2}(21) == 'N'
                SegmentedSampelsInfoLists.Task{i-2} = 'nback';
            else
                SegmentedSampelsInfoLists.Task{i-2} = 'ND';
            end
        else
            if contains(lower(SegmentedSampelsInfo{i-2}.Name),'bac')
                SegmentedSampelsInfoLists.Task{i-2} = 'nback';
            elseif contains(lower(SegmentedSampelsInfo{i-2}.Name),'odd')
                SegmentedSampelsInfoLists.Task{i-2} = 'oddball';
            else
                SegmentedSampelsInfoLists.Task{i-2} = 'both';
            end
        end
    end
end

end

function [SegmentedSampelsInfo, SegmentedSampelsInfoLists] = RawFilesInfo(Samplelisting)

SegmentedSampelsInfo = {};
for i = 3:length(Samplelisting)
    SegmentedSampelsInfo{i-2}.Name = Samplelisting(i).name;
    SegmentedSampelsInfo{i-2}.Address = [Samplelisting(i).folder '\' Samplelisting(i).name];
    if SegmentedSampelsInfo{i-2}.Name(1)== 'R'
        SegmentedSampelsInfo{i-2}.Name(1:2) = '';
    end
    switch SegmentedSampelsInfo{i-2}.Name(5)
        case 'F'
            SegmentedSampelsInfo{i-2}.Sport = 'football';
        case 'L'
            SegmentedSampelsInfo{i-2}.Sport = 'Lacrosse';
        case 'S'
            SegmentedSampelsInfo{i-2}.Sport = 'soccer';
        case 'R'
            SegmentedSampelsInfo{i-2}.Sport = 'rugby';
        otherwise
            SegmentedSampelsInfo{i-2}.Sport = 'ND';
    end
    
    if SegmentedSampelsInfo{i-2}.Name(1:3) == 'AFP'
        SegmentedSampelsInfo{i-2}.IDType = 'AFP';
        SegmentedSampelsInfo{i-2}.ID = SegmentedSampelsInfo{i-2}.Name(5:9);
        switch SegmentedSampelsInfo{i-2}.Name(10)
            case 'f'
                SegmentedSampelsInfo{i-2}.Gender = 'female';
            case 'm'
                SegmentedSampelsInfo{i-2}.Gender = 'male';
        end
        tmpy = split(SegmentedSampelsInfo{i-2}.Name,'_');
        SegmentedSampelsInfo{i-2}.Year = str2double(tmpy{3}(2:end));
        SegmentedSampelsInfo{i-2}.Session = str2double(cell2mat(strtrim(extractBetween(SegmentedSampelsInfo{i-2}.Name,'ses','_'))));
        if SegmentedSampelsInfo{i-2}.Name(16) == '_'
            if SegmentedSampelsInfo{i-2}.Name(17) == 'o' || SegmentedSampelsInfo{i-2}.Name(17) == 'O'
                SegmentedSampelsInfo{i-2}.Task = 'oddball';
            elseif SegmentedSampelsInfo{i-2}.Name(17) == 'n' || SegmentedSampelsInfo{i-2}.Name(17) == 'N'
                SegmentedSampelsInfo{i-2}.Task = 'nback';
            else
                SegmentedSampelsInfo{i-2}.Task = 'ND';
            end
        else
            if contains(lower(SegmentedSampelsInfo{i-2}.Name),'bac')
                SegmentedSampelsInfo{i-2}.Task = 'nback';
            elseif contains(lower(SegmentedSampelsInfo{i-2}.Name),'odd')
                SegmentedSampelsInfo{i-2}.Task = 'oddball';
            else
                SegmentedSampelsInfo{i-2}.Task = 'both';
            end
        end
    else
        SegmentedSampelsInfo{i-2}.IDType = 'TBIC';
        SegmentedSampelsInfo{i-2}.ID = SegmentedSampelsInfo{i-2}.Name(7:11);
        switch SegmentedSampelsInfo{i-2}.Name(12)
            case 'f'
                SegmentedSampelsInfo{i-2}.Gender = 'female';
            case 'm'
                SegmentedSampelsInfo{i-2}.Gender = 'male';
        end
        tmpy = split(SegmentedSampelsInfo{i-2}.Name,'_');
        SegmentedSampelsInfo{i-2}.Year = str2double(tmpy{3}(2:end));
        SegmentedSampelsInfo{i-2}.Session = str2double(cell2mat(strtrim(extractBetween(SegmentedSampelsInfo{i-2}.Name,'ses','_'))));
        if SegmentedSampelsInfo{i-2}.Name(21) == '_'
            if SegmentedSampelsInfo{i-2}.Name(22) == 'o' || SegmentedSampelsInfo{i-2}.Name(22) == 'O'
                SegmentedSampelsInfo{i-2}.Task = 'oddball';
            elseif SegmentedSampelsInfo{i-2}.Name(22) == 'n' || SegmentedSampelsInfo{i-2}.Name(22) == 'N'
                SegmentedSampelsInfo{i-2}.Task = 'nback';
            else
                SegmentedSampelsInfo{i-2}.Task = 'ND';
            end
        else
            if contains(lower(SegmentedSampelsInfo{i-2}.Name),'bac')
                SegmentedSampelsInfo{i-2}.Task = 'nback';
            elseif contains(lower(SegmentedSampelsInfo{i-2}.Name),'odd')
                SegmentedSampelsInfo{i-2}.Task = 'oddball';
            else
                SegmentedSampelsInfo{i-2}.Task = 'both';
            end
        end
    end
end

SegmentedSampelsInfoLists = {};
for i = 3:length(Samplelisting)
    SegmentedSampelsInfoLists.Name{i-2} = Samplelisting(i).name;
    SegmentedSampelsInfoLists.Address{i-2} = [Samplelisting(i).folder '\' Samplelisting(i).name];
    if SegmentedSampelsInfoLists.Name{i-2}(1)== 'R'
        SegmentedSampelsInfoLists.Name{i-2}(1:2) = '';
    end
    switch SegmentedSampelsInfoLists.Name{i-2}(5)
        case 'F'
            SegmentedSampelsInfoLists.Sport{i-2} = 'football';
        case 'L'
            SegmentedSampelsInfoLists.Sport{i-2} = 'Lacrosse';
        case 'S'
            SegmentedSampelsInfoLists.Sport{i-2} = 'soccer';
        case 'R'
            SegmentedSampelsInfoLists.Sport{i-2} = 'rugby';
        otherwise
            SegmentedSampelsInfoLists.Sport{i-2} = 'ND';
    end
    
    if SegmentedSampelsInfoLists.Name{i-2}(1:3) == 'AFP'
        SegmentedSampelsInfoLists.IDType{i-2} = 'AFP';
        SegmentedSampelsInfoLists.ID{i-2} = SegmentedSampelsInfoLists.Name{i-2}(5:9);
        switch SegmentedSampelsInfoLists.Name{i-2}(10)
            case 'f'
                SegmentedSampelsInfoLists.Gender{i-2} = 'female';
            case 'm'
                SegmentedSampelsInfoLists.Gender{i-2} = 'male';
        end
        tmpy = split(SegmentedSampelsInfoLists.Name{i-2},'_');
        SegmentedSampelsInfoLists.Year{i-2} = str2double(tmpy{3}(2:end));
        SegmentedSampelsInfoLists.Session{i-2} = str2double(cell2mat(strtrim(extractBetween(SegmentedSampelsInfoLists.Name{i-2},'ses','_'))));
        if SegmentedSampelsInfoLists.Name{i-2}(16) == '_'
            if SegmentedSampelsInfoLists.Name{i-2}(17) == 'o' || SegmentedSampelsInfoLists.Name{i-2}(17) == 'O'
                SegmentedSampelsInfoLists.Task{i-2} = 'oddball';
            elseif SegmentedSampelsInfoLists.Name{i-2}(17) == 'n' || SegmentedSampelsInfoLists.Name{i-2}(17) == 'N'
                SegmentedSampelsInfoLists.Task{i-2} = 'nback';
            else
                SegmentedSampelsInfoLists.Task{i-2} = 'ND';
            end
        else
            if contains(lower(SegmentedSampelsInfoLists.Name{i-2}),'bac')
                SegmentedSampelsInfoLists.Task{i-2} = 'nback';
            elseif contains(lower(SegmentedSampelsInfoLists.Name{i-2}),'odd')
                SegmentedSampelsInfoLists.Task{i-2} = 'oddball';
            else
                SegmentedSampelsInfoLists.Task{i-2} = 'both';
            end
        end
    else
        SegmentedSampelsInfoLists.IDType{i-2} = 'TBIC';
        SegmentedSampelsInfoLists.ID{i-2} = SegmentedSampelsInfoLists.Name{i-2}(7:11);
        switch SegmentedSampelsInfoLists.Name{i-2}(12)
            case 'f'
                SegmentedSampelsInfoLists.Gender{i-2} = 'female';
            case 'm'
                SegmentedSampelsInfoLists.Gender{i-2} = 'male';
        end
        tmpy = split(SegmentedSampelsInfoLists.Name{i-2},'_');
        SegmentedSampelsInfoLists.Year{i-2} = str2double(tmpy{3}(2:end));
        SegmentedSampelsInfoLists.Session{i-2} = str2double(cell2mat(strtrim(extractBetween(SegmentedSampelsInfoLists.Name{i-2},'ses','_'))));
        if SegmentedSampelsInfoLists.Name{i-2}(21) == '_'
            if SegmentedSampelsInfoLists.Name{i-2}(22) == 'o' || SegmentedSampelsInfoLists.Name{i-2}(22) == 'O'
                SegmentedSampelsInfoLists.Task{i-2} = 'oddball';
            elseif SegmentedSampelsInfoLists.Name{i-2}(22) == 'n' || SegmentedSampelsInfoLists.Name{i-2}(22) == 'N'
                SegmentedSampelsInfoLists.Task{i-2} = 'nback';
            else
                SegmentedSampelsInfoLists.Task{i-2} = 'ND';
            end
        else
            if contains(lower(SegmentedSampelsInfoLists.Name{i-2}),'bac')
                SegmentedSampelsInfoLists.Task{i-2} = 'nback';
            elseif contains(lower(SegmentedSampelsInfoLists.Name{i-2}),'odd')
                SegmentedSampelsInfoLists.Task{i-2} = 'oddball';
            else
                SegmentedSampelsInfoLists.Task{i-2} = 'both';
            end
        end
    end
end

end

function [EpochRaw, FreqBandRaw] = PartitioningEpochs(tmpPartitionName,Data, Fs, BlockSize, WindowMove, FreqBand, Format)

    tmpPartitionFile = ['D:\Concussion\Dataset\Partitions\' tmpPartitionName];

    if exist(tmpPartitionFile,'file')
        disp(['------------Loading ' tmpPartitionFile '  Partitions--------------------']);
        load(tmpPartitionFile)
    else
        disp('--------------Partitioning epochs ---------------------');
        [ EpochRaw, FreqBandRaw ] = PartinionData( Data, Fs, BlockSize, WindowMove, FreqBand, Format,'Tukey');
        disp(['------------Saving ' tmpPartitionFile '  Partitions--------------------']);
        save(tmpPartitionFile,'EpochRaw','FreqBandRaw');
    end

end

% * -1) Filters raw data 
function [ FilteredData ] = PreProccFilters( PureData,Fs,varargin )
%function [ FilteredData ] = HighPassFilter( PureData,Fs )
%   20th Order Cheby IIR filter BandPass filter in range of (0.2)Hz to 100 Hz
%   Also in this function the 60Hz (59.5 to 60.5) noise is filtered by an Elliptic filter

 %%---Filtering Data
order    = 20;
fpass  = (0.2); % default
if ~isempty(varargin) && isa(varargin{1},'double')
    fpass  = varargin{1};
end

% [b,a] = butter(order,fpass, 'high');			   

% hpFilt = designfilt('highpassiir', 'FilterOrder', order, 'PassbandFrequency', 0.1, ...
%                     'StopbandAttenuation', 40, 'PassbandRipple', 0.05, ...
%                     'SampleRate', Fs, 'DesignMethod', 'ellip');
		
% lpFilt = designfilt('lowpassfir','PassbandFrequency',0.8, ...
%                     'StopbandFrequency',0.85,'PassbandRipple',0.1, ...
%                     'StopbandAttenuation',65,'DesignMethod','kaiserwin');
                
lpFilt = designfilt('lowpassfir','PassbandFrequency',0.64, ...
                    'StopbandFrequency',0.7,'PassbandRipple',0.1, ...
                    'StopbandAttenuation',65,'DesignMethod','kaiserwin');
                 
% NotchFilt = designfilt('bandstopiir', 'FilterOrder', order, 'StopbandAttenuation', 40, ...
%                'StopbandFrequency1', 59.5, 'StopbandFrequency2', 60.5,  ...
%                'SampleRate', Fs, 'DesignMethod', 'cheby2');
           
hpFilt = designfilt('highpassfir','StopbandFrequency',0.004, ...
                    'PassbandFrequency',0.008,'PassbandRipple',0.5, ...
                    'StopbandAttenuation',65,'DesignMethod','kaiserwin');
                
bsFilt = designfilt('bandstopfir','FilterOrder',250, ...
                    'CutoffFrequency1',59.5,'CutoffFrequency2',60.5, ...
                    'SampleRate',Fs);
FilteredData = zeros(size(PureData,1),size(PureData,2));
for i=1:size(PureData,1)
    tmp = PureData(i,:);
    tmpbs = filtfilt(bsFilt,tmp);
    tmphp = filtfilt(hpFilt,tmpbs);
    FilteredData(i,:) = filtfilt(lpFilt,tmphp);
end
end

% * 0) Given the Mastersheet the => samples info are generating
function [Sample,T] = MasterSheet2Sample(MasterSheetAdd,SheetName,RawSampels,FilteredSampels,CleanSampels)
    [~,~,raw] = xlsread(MasterSheetAdd,SheetName);
    raw(:,[8 10 15:18])=[];

    T = table(raw(2:end,1),raw(2:end,2),raw(2:end,3),raw(2:end,4),raw(2:end,5),raw(2:end,6),...
              raw(2:end,7),raw(2:end,8),raw(2:end,9),raw(2:end,10),raw(2:end,11),raw(2:end,12),raw(2:end,13),...
              'VariableNames',{'ID','AFPID','YearOfStudy','Session','SessionComment','Sport','Gender','Conc_hx',...
                               'MidSession_Session','DateOfTest','DateOfBirth','Age','BadElectrodes'});

    for t = 1:size(T,1)
        if isnan(cell2mat(T{t,1}))
            T(t:end,:) = [];
            break;
        end
    end
    
    Sample = cell(1,size(T,1));
    for s = 1:size(T,1)      
        tmpsample = table2struct(T(s,:));
        tmpsample.RawFiles = {};
        for i = 1:length(RawSampels)
            if strcmp(RawSampels{i}.IDType,'TBIC')
                if str2num(RawSampels{i}.ID) == T.ID{s}
                    if RawSampels{i}.Session == T.Session{s}
                        if RawSampels{i}.Year == str2double(T.YearOfStudy{s}(1))
                            tmpsample.RawFiles{end+1}.Task = RawSampels{i}.Task;
                            tmpsample.RawFiles{end}.Address = RawSampels{i}.Address;
                        end
                    end
                end
            end
            if ~isnan(T.AFPID{s}) && strcmp(RawSampels{i}.IDType,'AFP')
                if str2num(RawSampels{i}.ID) == T.AFPID{s}
                    if RawSampels{i}.Session == T.Session{s}
                        if RawSampels{i}.Year == str2double(T.YearOfStudy{s}(1))
                            tmpsample.RawFiles{end+1}.Task = RawSampels{i}.Task;
                            tmpsample.RawFiles{end}.Address = RawSampels{i}.Address;
                        end
                    end
                end               
            end     
        end
        tmpsample.FilteredFiles = {};
        for i = 1:length(FilteredSampels)
            if strcmp(FilteredSampels{i}.IDType,'TBIC')
                if str2num(FilteredSampels{i}.ID) == T.ID{s}
                    if FilteredSampels{i}.Session == T.Session{s}
                        if FilteredSampels{i}.Year == str2double(T.YearOfStudy{s}(1))
                            tmpsample.FilteredFiles{end+1}.Task = FilteredSampels{i}.Task;
                            tmpsample.FilteredFiles{end}.Address = FilteredSampels{i}.Address;
                        end
                    end
                end
            end
            if ~isnan(T.AFPID{s}) && strcmp(FilteredSampels{i}.IDType,'AFP')
                if str2num(FilteredSampels{i}.ID) == T.AFPID{s}
                    if FilteredSampels{i}.Session == T.Session{s}
                        if FilteredSampels{i}.Year == str2double(T.YearOfStudy{s}(1))
                            tmpsample.FilteredFiles{end+1}.Task = FilteredSampels{i}.Task;
                            tmpsample.FilteredFiles{end}.Address = FilteredSampels{i}.Address;
                        end
                    end
                end               
            end     
        end
        tmpsample.CleanFiles = {};        
        for i = 1:length(CleanSampels)
            if strcmp(CleanSampels{i}.IDType,'TBIC')
                if str2num(CleanSampels{i}.ID) == T.ID{s}
                    if CleanSampels{i}.Session == T.Session{s}
                        if CleanSampels{i}.Year == str2double(T.YearOfStudy{s}(1))
                            tmpsample.CleanFiles{end+1}.Task = CleanSampels{i}.Task;
                            tmpsample.CleanFiles{end}.Address = CleanSampels{i}.Address;
                        end
                    end
                end
            end
            if ~isnan(T.AFPID{s}) && strcmp(CleanSampels{i}.IDType,'AFP')
                if str2num(CleanSampels{i}.ID) == T.AFPID{s}
                    if CleanSampels{i}.Session == T.Session{s}
                        if CleanSampels{i}.Year == str2double(T.YearOfStudy{s}(1))
                            tmpsample.CleanFiles{end+1}.Task = CleanSampels{i}.Task;
                            tmpsample.CleanFiles{end}.Address = CleanSampels{i}.Address;
                        end
                    end
                end               
            end     
        end
        Sample{s} = tmpsample;
    end
end


% * 1) FreqBand initials
function [cCenters] = ElecMap()

Centers = [1044.37878787879,279.900000000000;975.106060606061,278.360606060606;910.451515151515,287.596969696970;862.730303030303,315.306060606060;818.087878787879,343.015151515151;791.918181818182,386.118181818182;759.590909090909,427.681818181818;728.803030303030,470.784848484848;704.172727272727,523.124242424242;962.790909090909,210.627272727273;898.136363636364,226.021212121212;842.718181818182,252.190909090909;799.615151515152,287.596969696970;756.512121212121,330.700000000000;730.342424242424,375.342424242424;698.015151515152,426.142424242424;676.463636363636,484.639393939394;879.663636363636,165.984848484848;821.166666666667,198.312121212121;770.366666666667,235.257575757576;728.803030303030,284.518181818182;701.093939393939,330.700000000000;671.845454545455,383.039393939394;650.293939393940,444.615151515151;782.681818181818,144.433333333333;730.342424242424,189.075757575758;690.318181818182,235.257575757576;659.530303030303,289.136363636364;642.596969696970,343.015151515151;622.584848484849,412.287878787879;730.342424242424,95.1727272727271;676.463636363636,144.433333333333;642.596969696970,198.312121212121;616.427272727273,250.651515151515;599.493939393940,315.306060606060;599.493939393940,373.803030303030;585.639393939394,167.524242424242;564.087878787879,229.100000000000;547.154545454545,287.596969696970;537.918181818182,358.409090909091;553.312121212121,416.906060606060;585.639393939394,452.312121212121;613.348484848485,490.796969696970;647.215151515151,532.360606060606;687.239393939394,575.463636363636;499.433333333333,209.087878787879;484.039393939394,279.900000000000;482.500000000000,343.015151515151;493.275757575758,412.287878787879;516.366666666667,469.245454545454;553.312121212121,503.112121212121;594.875757575758,544.675757575757;644.136363636364,587.778787878788;420.924242424242,282.978787878788;419.384848484849,358.409090909091;434.778787878788,423.063636363636;462.487878787879,486.178787878788;504.051515151515,523.124242424242;548.693939393939,561.609090909091;596.415151515152,604.712121212121;354.730303030303,383.039393939394;373.203030303030,453.851515151515;413.227272727273,507.730303030303;450.172727272727,553.912121212121;500.972727272727,581.621212121212;557.930303030303,620.106060606060;282.378787878788,418.445454545454;322.403030303030,487.718181818182;357.809090909091,555.451515151515;405.530303030303,592.396969696970;456.330303030303,620.106060606060;510.209090909091,635.500000000000;223.881818181818,469.245454545454;367.045454545455,644.736363636364;424.003030303030,673.984848484849;482.500000000000,687.839393939394;536.378787878788,687.839393939394;581.021212121212,667.827272727273;627.203030303030,652.433333333333;671.845454545455,629.342424242424;730.342424242424,607.790909090909;173.081818181818,555.451515151515;339.336363636364,704.772727272727;397.833333333333,729.403030303030;460.948484848485,743.257575757576;514.827272727273,737.100000000000;568.706060606061,732.481818181818;622.584848484849,709.390909090909;670.306060606061,687.839393939394;728.803030303030,650.893939393939;128.439393939394,629.342424242424;186.936363636364,687.839393939394;248.512121212121,737.100000000000;308.548484848485,757.112121212121;380.900000000000,789.439393939394;447.093939393940,801.754545454546;505.590909090909,798.675757575757;564.087878787879,784.821212121212;613.348484848485,763.269696969697;674.924242424242,737.100000000000;733.421212121212,743.257575757576;163.845454545455,777.124242424242;234.657575757576,812.530303030303;297.772727272727,841.778787878788;371.663636363636,860.251515151515;445.554545454545,869.487878787879;514.827272727273,864.869696969697;573.324242424243,849.475757575758;631.821212121212,821.766666666667;670.306060606061,787.900000000000;226.960606060606,897.196969696970;303.930303030303,921.827272727273;377.821212121212,940.300000000000;457.869696969697,937.221212121212;527.142424242424,929.524242424242;588.718181818182,907.972727272727;647.215151515151,880.263636363636;696.475757575758,838.700000000000;730.342424242424,794.057575757576;317.784848484849,1006.49393939394;399.372727272727,1009.57272727273;477.881818181818,1008.03333333333;554.851515151515,995.718181818182;624.124242424242,966.469696969697;682.621212121212,932.603030303030;730.342424242424,889.500000000000;761.130303030303,838.700000000000;787.300000000000,787.900000000000;785.760606060606,737.100000000000;787.300000000000,687.839393939394;785.760606060606,633.960606060606;773.445454545455,572.384848484848;430.160606060606,1089.62121212121;519.445454545455,1077.30606060606;596.415151515152,1060.37272727273;664.148484848485,1024.96666666667;730.342424242424,980.324242424242;776.524242424242,929.524242424242;813.469696969697,878.724242424242;831.942424242424,823.306060606060;845.796969696970,761.730303030303;839.639393939394,710.930303030303;836.560606060606,652.433333333333;815.009090909091,587.778787878788;579.481818181818,1146.57878787879;653.372727272727,1114.25151515152;730.342424242424,1077.30606060606;793.457575757576,1024.96666666667;839.639393939394,969.548484848485;871.966666666667,911.051515151515;887.360606060606,847.936363636364;895.057575757576,784.821212121212;893.518181818182,730.942424242424;878.124242424242,667.827272727273;861.190909090909,604.712121212121;807.312121212121,1109.63333333333;862.730303030303,1060.37272727273;905.833333333333,992.639393939394;932.003030303030,929.524242424242;945.857575757576,864.869696969697;958.172727272727,794.057575757576;948.936363636364,740.178787878788;925.845454545455,689.378787878788;905.833333333333,624.724242424242;887.360606060606,1145.03939393939;941.239393939394,1077.30606060606;981.263636363636,1008.03333333333;1005.89393939394,938.760606060606;1016.66969696970,869.487878787879;1015.13030303030,801.754545454546;1004.35454545455,741.718181818182;981.263636363636,687.839393939394;952.015151515152,638.578787878788;1032.06363636364,1088.08181818182;1062.85151515152,1009.57272727273;1081.32424242424,937.221212121212;1087.48181818182,858.712121212121;1081.32424242424,794.057575757576;1062.85151515152,730.942424242424;1038.22121212121,673.984848484849;1005.89393939394,618.566666666667;962.790909090909,581.621212121212;915.069696969697,561.609090909091;867.348484848485,544.675757575757;815.009090909091,529.281818181818;751.893939393940,523.124242424242;1142.90000000000,1001.87575757576;1156.75454545455,918.748484848485;1159.83333333333,843.318181818182;1153.67575757576,757.112121212121;1122.88787878788,703.233333333333;1092.10000000000,644.736363636364;1055.15454545455,590.857575757576;1010.51212121212,550.833333333333;959.712121212121,523.124242424242;904.293939393940,503.112121212121;845.796969696970,484.639393939394;782.681818181818,484.639393939394;1236.80303030303,892.578787878788;1226.02727272727,812.530303030303;1213.71212121212,737.100000000000;1102.87575757576,553.912121212121;1048.99696969697,507.730303030303;999.736363636364,480.021212121212;942.778787878788,463.087878787879;873.506060606061,447.693939393939;808.851515151515,447.693939393939;1298.37878787879,774.045454545455;1273.74848484849,684.760606060606;1144.43939393939,487.718181818182;1085.94242424242,449.233333333333;1027.44545454545,421.524242424242;967.409090909091,412.287878787879;910.451515151515,412.287878787879;835.021212121212,410.748484848485;1335.32424242424,624.724242424242;1289.14242424242,553.912121212121;1233.72424242424,460.009090909091;1178.30606060606,415.366666666667;1107.49393939394,379.960606060606;1042.83939393939,353.790909090909;978.184848484848,347.633333333333;919.687878787879,356.869696969697;862.730303030303,373.803030303030;1167.53030303030,353.790909090909;1139.82121212121,293.754545454545;1229.10606060606,367.645454545454;1295.30000000000,426.142424242424;1356.87575757576,486.178787878788;1105.95454545455,218.324242424242;1207.55454545455,269.124242424242;1287.60303030303,315.306060606060;1353.79696969697,373.803030303030;1050.53636363636,138.275757575757;1153.67575757576,156.748484848485;1244.50000000000,196.772727272727;1321.46969696970,247.572727272727;1001.27575757576,81.3181818181815;1101.33636363636,76.6999999999998;1186.00303030303,96.7121212121210;457.869696969697,81.3181818181815;360.887878787879,81.3181818181815;274.681818181818,101.330303030303;410.148484848485,141.354545454545;307.009090909091,161.366666666667;217.724242424242,198.312121212121;136.136363636364,252.190909090909;354.730303030303,222.942424242424;253.130303030303,273.742424242424;173.081818181818,319.924242424242;106.887878787879,378.421212121212;319.324242424242,298.372727272727;293.154545454545,358.409090909091;231.578787878788,375.342424242424;163.845454545455,430.760606060606;105.348484848485,490.796969696970];


cXs = (-1)* (Centers(:,1) - mean(Centers(:,1)));
cYs = (-1)* (Centers(:,2) - mean(Centers(:,2)));
cCenters = [cXs cYs];   % Corrected Centers, transfer the origin to ~(0,0)

end

function [Data, Electrodes, TestPeriods] = QualityControl(TMPrecorddata)

    Electrodes = [1:size(TMPrecorddata,1)];
        
    tmpZeros = find(sum(TMPrecorddata)==0);
    tmpStart = find(diff(tmpZeros)>1);
    
    if length(tmpStart) == 1
        TestPeriods{1} = [tmpZeros(tmpStart)+1 tmpZeros(tmpStart+1)-1];
        tmpData = TMPrecorddata(:,TestPeriods{1}(1):TestPeriods{1}(2));
    elseif isempty(tmpStart)
        if isempty(tmpZeros)
            TestPeriods{1} = [1 length(TMPrecorddata)];
            tmpData = TMPrecorddata(:,TestPeriods{1}(1):TestPeriods{1}(2));
        else
            if length([1:tmpZeros(1)-1]) > length([tmpZeros(end)+1:length(TMPrecorddata)])
                TestPeriods{1} = [1 tmpZeros(1)-1];
                tmpData = TMPrecorddata(:,TestPeriods{1}(1):TestPeriods{1}(2));
            else
                TestPeriods{1} = [tmpZeros(end)+1 length(TMPrecorddata)];
                tmpData = TMPrecorddata(:,TestPeriods{1}(1):TestPeriods{1}(2));
            end
        end
    else
        TestPeriods{1} = [tmpZeros(tmpStart(2))+1 tmpZeros(tmpStart(2)+1)-1];
        tmpData = TMPrecorddata(:,TestPeriods{1}(1):TestPeriods{1}(2));
    end

    
    SumFalses = sum((tmpData == -5000),2) + sum((tmpData == 5000),2);
    tmpBrokeChann1 = find(SumFalses>20);

    % Consecutive Disconnection
    
    tmpDiff = diff(tmpData,1,2);
    tmpBrokeChann2 = [];
    for j = 1:size(tmpData,1)
        if ismember(j,tmpBrokeChann1)
            continue;
        end
        maxZerro = 1;
        currentRun = 1;
        for i = 1:numel(tmpDiff(j,:))
            if tmpDiff(j,i)==0
                currentRun = currentRun + 1;
                maxZerro = max(maxZerro,currentRun);
            else
                currentRun = 1;
            end
        end
        if maxZerro > 250
            tmpBrokeChann2 = [tmpBrokeChann2 j];
        end
    end

%     tmpBrokeChann2 = find(sum(tmpData(:,1:end-1)-tmpData(:,2:end)==0,2) > 250);  
    
    BrokenChann = union(tmpBrokeChann1,tmpBrokeChann2);
    
    Electrodes(BrokenChann) = [];
    
    disp(['Broken Chann1 (Completely disconnected) = ' num2str(length(tmpBrokeChann1)) ' Broken Chann2 (disconnected at least for a second) = ' num2str(length(tmpBrokeChann2))]);
    
    Data = tmpData(Electrodes,:);

    % figure;
    % plot((tmpData-mean(tmpData,2))')
    
    
    % figure;
    % plot((tmpData(SubZeros))')

%     tmpData = [];
%     tmpData = recorddata(Electrodes,:);

%     figure;
%     hist(SumFalses,100)
%     title('Broken Chann. Step 1')


end
% Functions inside this code
% * 1) FreqBand initials
function FreqBand = FreqBand_init(varargin)
    %- MaxF is the Nyquist freq.
    
     if ~isempty(varargin) && isa(varargin{1},'double')
         MaxF = floor(varargin{1});
     else
         MaxF = 125;
     end
    
    FreqBand(1).Range = [0.5 4];
    FreqBand(1).Name = 'Delta';
    FreqBand(2).Range = [4 8];
    FreqBand(2).Name = 'Theta';
    FreqBand(3).Range = [8 14];
    FreqBand(3).Name = 'Alpha';
    FreqBand(4).Range = [14 30];
    FreqBand(4).Name = 'Beta';
    FreqBand(5).Range = [30 80];
    FreqBand(5).Name = 'Gamma';
    FreqBand(6).Range = [80 MaxF];
    FreqBand(6).Name = 'High-Gamma';
    FreqBand(7).Range = [0.5 8];
    FreqBand(7).Name = 'Delta-Theta'; 
    FreqBand(8).Range = [4 14];
    FreqBand(8).Name = 'Theta-Alpha'; 
    FreqBand(9).Range = [8 30];
    FreqBand(9).Name = 'Alpha-Beta';
    FreqBand(10).Range = [14 80];
    FreqBand(10).Name = 'Beta-Gamma';
    FreqBand(11).Range = [30 MaxF];
    FreqBand(11).Name = 'Gammas';
    FreqBand(12).Range = [0.5 MaxF];
    FreqBand(12).Name = 'All';
    
end

% * 2) Reading EEG file in ".edf" format
function [Data, Fs, Labels, recorddata, header] = ReadEDFformatEEG(FileAddress, Electrodes)
    [header, recorddata] = edfread(FileAddress);
    Fs=(header.samples(1))/(header.duration);
    %%---Labels
    for i=1:length(Electrodes)
        Labels{i} = header.label{Electrodes(i)};
    end
    Data = HighPassFilter( recorddata,Fs);
end

% * 3) Reading EEG file in ".xlsx" format
function [Data, Fs, Labels] = ReadExcelformatEEG(FileAddress)
    Fs = 256;  
    [num,txt,~] = xlsread(FileAddress);
    %%---Labels
    Labels = txt;
    Data = HighPassFilter( num',Fs );
end

% * 4) Initiate the connection network matrices (Zero matrices)
function [FreqBand] = Connection_init(FreqBand, RegionLabels)
    for f = 1:length(FreqBand)
        FreqBand(f).VectorConn = zeros(length(FreqBand(f).Energies),length(RegionLabels));
        FreqBand(f).VectorConnstr = zeros(length(FreqBand(f).Energies),length(RegionLabels));
        FreqBand(f).VectorConnDist = zeros(length(FreqBand(f).Energies),length(RegionLabels));
        FreqBand(f).VectorConnstrDist = zeros(length(FreqBand(f).Energies),length(RegionLabels));
    end
end

% * 5) Calculating the Normalized epoch Total energies of an epoch
%      base on the [2N to N] previous epochs 
function [FreqBand] = NormalizedEnergies(FreqBand, f, i, N)

% FreqBand(f).TotEnergiesNorm(i) =FreqBand(f).TotEnergies(i)/mean(FreqBand(f).TotEnergies(i-2*N:i-N));
% FreqBand(f).TotEnergiesNorm(i) =((FreqBand(f).TotEnergies(i)/mean(FreqBand(f).TotEnergies(i-2*N:i-N)))/mean(FreqBand(f).TotEnergiesNorm(1:2*N)));
% FreqBand(f).TotEnergiesNorm(i) =((FreqBand(f).TotEnergies(i)/mean(FreqBand(f).TotEnergies(i-2*N:i-N)))- min(FreqBand(f).TotEnergiesNorm(max(1,i-N):i-1)))/(max(FreqBand(f).TotEnergiesNorm(max(1,i-N):i-1))-min(FreqBand(f).TotEnergiesNorm(max(1,i-N):i-1)));
% FreqBand(f).TotEnergiesNorm(i) =((FreqBand(f).TotEnergies(i)/mean(FreqBand(f).TotEnergies(i-2*N:i-N)))- min(FreqBand(f).TotEnergiesNorm(max(1,i-2*N):i-1)))/(max(FreqBand(f).TotEnergiesNorm(max(1,i-2*N):i-1))-min(FreqBand(f).TotEnergiesNorm(max(1,i-2*N):i-1)));
% FreqBand(f).TotEnergiesNorm(i) =((FreqBand(f).TotEnergies(i)/mean(FreqBand(f).TotEnergies(i-2*N:i-N)))- min(FreqBand(f).TotEnergies(max(1,i-2*N):i-1)))/(max(FreqBand(f).TotEnergies(max(1,i-2*N):i-1))-min(FreqBand(f).TotEnergies(max(1,i-2*N):i-1)));
FreqBand(f).TotEnergiesNorm(i) =(FreqBand(f).TotEnergies(i)- min(FreqBand(f).TotEnergies(max(1,i-2*N):i-N)))/(max(FreqBand(f).TotEnergies(max(1,i-2*N):i-N))-min(FreqBand(f).TotEnergies(max(1,i-2*N):i-N)));


end

% * 6) Calculating the Average Total Energies of epochs 
%      base on the half of electrodes with less energies in each Freq band
function [FreqBand] = TotalAvgSqrEngery(FreqBand)
%          FreqBand(f).TotEnergies = mean(FreqBand(f).Energies.^2,2);
    for f = 1:length(FreqBand)
        tmp = sort(FreqBand(f).Energies.^2,2);
        FreqBand(f).TotEnergies = mean(tmp(:,1:ceil(end/2)),2);  
    end
end
function [FreqBand] = TotalAvgSqrEngery1(FreqBand)
%          FreqBand(f).TotEnergies = mean(FreqBand(f).Energies.^2,2);
    for f = 1:length(FreqBand)
        tmp = sort(FreqBand(f).Energies,2);
        FreqBand(f).TotEnergies = mean(tmp(:,1:ceil(end/(4/3))),2);  
    end
end

% * 7) Picking the FFT or DCT energies and the Filtered bands 
%      base on the defined "Format"
function [FreqBand] = PickVariables(FreqBand, Format)
    if Format == 1 || Format == 3
        for f = 1:length(FreqBand)
            FreqBand(f).Energies = FreqBand(f).FFTEnergies;
            FreqBand(f).FilteredDATA = FreqBand(f).FilteredDataFFT;
        end
    elseif Fromat == 2
        for f = 1:length(FreqBand)
            FreqBand(f).Energies = FreqBand(f).DCTEnergies;
            FreqBand(f).FilteredDATA = FreqBand(f).FilteredDataDCT;
        end
    end
end

%%
