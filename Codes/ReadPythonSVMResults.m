%% Concussion Analysis - Gathering Python SVM results
close all;
clear;
clc;

addpath(genpath('..\Functions'));
addpath(genpath('..\Functions\FastICA_25'));

%%
DistMtric = {'Euclid','Cos','Corr','AllFeat','AllFeatSegNorm','AllwPo','AllwPoSegNorm','Po2Euclid2','Po2Euclid3','Po3Euclid2','Po2Cos2','Po2Cos3','Po3Cos2','Po2Corr2','Po2Corr3','Po3Corr2'};
Regions = Region_Init_TBI();

FreqBand = FreqBand_init;
FreqBand([5,6,10,11]) = [];

PtOrder = [10004,10034,10045,10055,10073,10092,10048,10051,10058,10059,10087,10095];

RegionsName = {Regions.Name};

%%
XLLineAll = 1;
XLNameAll = ['Results\Results-TBI-Dist-Regions-SVM-',datestr(now,'HH-MM-SS'),'-Py.xlsx'];

% FolderAdd = 'D:\Concussion\07252020\Results\SVM results\';
% FolderAdd = 'I:\hard drive backup august 11-2020\Concussion\07312020\PythonSVMResults\';
% FolderAdd = 'D:\Concussion\09052020\PythonSVMResults\PythonSVMResults-HighRankedSegs\';
% FolderAdd = 'D:\Concussion\09052020\PythonSVMResults\PythonSVMResults-TimeRankedSegs\';
% FolderAdd = 'D:\Concussion\09052020\PythonSVMResults\PythonSVMResults-FreqRankedSegs\';
% FolderAdd = 'D:\Concussion\10082020\SVMResultsWAvgRmvWBadSegsRmv\Results_FreqRanked\';
% FolderAdd = 'D:\Concussion\10082020\SVMResultsWAvgRmvWBadSegsRmv\Results_TimeRanked\';
% FolderAdd = 'D:\Concussion\10082020\SVMResultsWAvgRmvWBadSegsRmv\Results_TimeFreqRanked\';
%FolderAdd = 'D:\Concussion\10082020\SVMResultsWAvgRmvWBadSegsRmv\Results_UnRanked\';

FolderAdd = 'D:\Concussion\10082020\All12PtTime04\Results_TimeRanked\';

for mett = 1:length(DistMtric)
    

    tmpOut = {'','','','','','','','','','','','','','','','','';
              'Freq Band','Region','ID10004','ID10034','ID10045','ID10055','ID10073','ID10092','ID10048','ID10051','ID10058','ID10059','ID10087','ID10095','Avg Concussed','Avg Control','Total Avg'};

    for f = 1:length(FreqBand)

        tmpAllACC = zeros(length(Regions),length(PtOrder));
        for pt = 1:length(PtOrder)
            
            tmptxtFileAdd = ['SVM_Dist_Freq_' DistMtric{mett} '_' FreqBand(f).Name '_ID_' num2str(PtOrder(pt)) '.txt'];

            tmpfile = tdfread([FolderAdd tmptxtFileAdd],'\t');
            
            tmpreg = cellstr([tmpfile.Regions(1:length(Regions),:)]);
            [~,~,idxsIntoB] = intersect(RegionsName,tmpreg,'stable');

            tmpSortedACC = tmpfile.ACC(idxsIntoB);
            
            tmpAllACC(:,pt) = tmpSortedACC;
        end

        tmpAllACC(:,end+1) = mean(tmpAllACC(:,1:6),2);
        tmpAllACC(:,end+1) = mean(tmpAllACC(:,7:12),2);
        tmpAllACC(:,end+1) = mean(tmpAllACC(:,1:12),2);

        tmpAllACCcell = num2cell(tmpAllACC);
        tmpres = {};
        
        for r = 1:length(Regions)
            tmpres(r,:) = {FreqBand(f).Name,RegionsName{r},tmpAllACCcell{r,:}};
        end
        
        tmpOut = [tmpOut;tmpres];

    end
    xlswrite(XLNameAll,tmpOut,DistMtric{mett},'A1');

end


%%
disp('-----Done!!-----');


%% Functions

function [Regions] = Region_Init_TBI()
    Regions = [];
    Regions(end+1).Name = 'R-Frontal';
    Regions(end).Indx = [1,2,3,4,5,6,10,11,12,13,14,15,18,19,20,21,25,26,222,223,224];
    Regions(end+1).Name = 'L-Frontal';
    Regions(end).Indx = [21,22,23,26,27,28,29,32,33,34,35,36,37,38,39,40,46,47,48,54];
    Regions(end+1).Name = 'Frontal';
    Regions(end).Indx = [1,2,3,4,5,6,10,11,12,13,14,15,18,19,20,21,25,26,222,223,224,22,23,26,27,28,29,32,33,34,35,36,37,38,39,40,46,47,48,54];
    Regions(end+1).Name = 'R-Temporal';
    Regions(end).Indx = [170,171,172,178,179,180,181,190,191,192,193,194,202,203,204,210,211,212,220,221];
    Regions(end+1).Name = 'L-Temporal';
    Regions(end).Indx = [55,56,57,61,62,63,64,68,69,70,71,74,75,76,83,84,85,94,95,96];
    Regions(end+1).Name = 'Temporal';
    Regions(end).Indx = [170,171,172,178,179,180,181,190,191,192,193,194,202,203,204,210,211,212,220,221,55,56,57,61,62,63,64,68,69,70,71,74,75,76,83,84,85,94,95,96];
    Regions(end+1).Name = 'R-Central';
    Regions(end).Indx = [7,8,81,90,130,131,132,142,143,144,154,155,163,164,173,182,183,184,185,186,195,196,197,198,205,206,207,213,214,215];
    Regions(end+1).Name = 'L-Central';
    Regions(end).Indx = [8,9,16,17,24,30,41,42,43,44,45,49,50,51,52,53,58,59,60,65,66,72,77,78,79,80,81,88,89,90];
    Regions(end+1).Name = 'Central';
    Regions(end).Indx = [7,8,81,90,130,131,132,142,143,144,154,155,163,164,173,182,183,184,185,186,195,196,197,198,205,206,207,213,214,215,9,16,17,24,30,41,42,43,44,45,49,50,51,52,53,58,59,60,65,66,72,77,78,79,80,88,89];
    Regions(end+1).Name = 'Parietal';
    Regions(end).Indx = [85,86,87,97,98,99,100,101,108,109,110,118,119,127,128,129,140,141,151,152,153,161,162,171];
    Regions(end+1).Name = 'Occipital';
    Regions(end).Indx = [106,107,108,114,115,116,117,123,124,125,126,136,137,138,139,147,148,149,150,151,158,159,160,168,169];
    Regions(end+1).Name = 'Right';
    Regions(end).Indx = [1,2,3,4,5,6,7,8,10,11,12,13,14,15,18,19,20,21,25,26,81,90,101,119,126,127,128,129,130,131,132,137,138,139,140,141,142,143,144,148,149,150,151,152,153,154,155,158,159,160,161,162,163,164,168,169,170,171,172,173,178,179,180,181,182,183,184,185,186,190,191,192,193,194,195,196,197,198,202,203,204,205,206,207,210,211,212,213,214,215,220,221,222,223,224];
    Regions(end+1).Name = 'Left';
    Regions(end).Indx = [8,9,16,17,21,22,23,24,26,27,28,29,30,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,68,69,70,71,72,74,75,76,77,78,79,80,81,83,84,85,86,87,88,89,90,94,95,96,97,98,99,100,101,106,107,108,109,110,114,115,116,117,118,123,124,125,126,136,137];
    Regions(end+1).Name = 'All';
    Regions(end).Indx = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,68,69,70,71,72,74,75,76,77,78,79,80,81,83,84,85,86,87,88,89,90,94,95,96,97,98,99,100,101,105,106,107,108,109,110,114,115,116,117,118,119,123,124,125,126,127,128,129,130,131,132,136,137,138,139,140,141,142,143,144,147,148,149,150,151,152,153,154,155,158,159,160,161,162,163,164,168,169,170,171,172,173,177, 178,179,180,181,182,183,184,185,186,190,191,192,193,194,195,196,197,198,202,203,204,205,206,207,210,211,212,213,214,215,220,221,222,223,224];
    Regions(end+1).Name = 'All 256 ch.';
    Regions(end).Indx = [1:256];
end

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
