% DOT-HUB toolbox Display data.
%
% Georgina's processing script, v0.0001;

%Loop around subjects (.LUMO directories);
%cfgFileName = [filepath '/ExampleData/Example1/preproPipelineExample1.cfg'];

fileList = dir('*.LUMO');

plotFlag = 0;

for i = 1:length(fileList)
    LUMODirName = fileList(i).name;
    
    % Covert .LUMO to .nirs
    [nirs, nirsFileName, SD3DFileName] = DOTHUB_LUMO2nirs(LUMODirName);
    
    % Run data quality checks - this produces multiple figures, so comment out for speed.
    if plotFlag
        DOTHUB_dataQualityCheck(nirsFileName);
    end
    %(include printFigFlag == 1 to output all the figures)
    %disp('Examine data quality figures, press any key to continue');
    %pause
    
    % Run Homer2 pre-processing pipeline using .cfg file. Alternatively you can run line by line (as per commented below).
    %[prepro, preproFileName] = DOTHUB_runHomerPrepro(nirsFileName,cfgFileName);
    
    %%%%Equivalent line-by-line Homer2 calls and prepro write:
    dod = hmrIntensity2OD(nirs.d);
    SD3D = enPruneChannels(nirs.d,nirs.SD3D,ones(size(nirs.t)),[0 1e11],12,[0 100],0);
    
    %Force MeasListAct to be the same across wavelengths
    SD3D = DOTHUB_balanceMeasListAct(SD3D);
    
    %Set SD2D
    SD2D = nirs.SD;
    SD2D.MeasListAct = SD3D.MeasListAct;
    
    %Convert to DOD
    dod = hmrBandpassFilt(dod,nirs.t,0.01,0.25);
    dc = hmrOD2Conc(dod,SD3D,[6 6]);
    dc = dc*1e6; %Homer works in Molar by default, we use uMolar.
    
    %Block avg (no regression)
    [dcAvg,dcAvgStd,tHRF] = hmrBlockAvg(dc,nirs.s,nirs.t,[-15 45]);
 
    %Regress short channels
    %Full time course regression. Overfit?
    
    dc_reg = DOTHUB_hmrSSRegressionByChannel(dc,SD3D,12,1); %This is a custom SS regression script.
    
    %Deconv model
    %[dcAvg_reg,dcAvgStd_reg,~] = hmrBlockAvg(dc_reg,nirs.s,nirs.t,[-10 50]);
    
    [dcAvg_ded,dcAvgStd_ded,~] = hmrDeconvHRF_DriftSS(dc,nirs.s,nirs.t,SD3D,[],ones(size(nirs.t)),[-15 45],1,1,[1 1],12,1,0,0);
    
    %Block avg (with regression)
    [dcAvg_reg,dcAvgStd_reg,~] = hmrBlockAvg(dc_reg,nirs.s,nirs.t,[-15 45]);

    %Output the data for every run of loop
    dcAvg_group(i,:,:,:) = dcAvg;
    dcAvgStd_group(i,:,:,:) = dcAvgStd;
    
    dcAvg_reg_group(i,:,:,:) = dcAvg_reg;
    dcAvgStd_reg_group(i,:,:,:) = dcAvgStd_reg;
    
    dcAvg_ded_group(i,:,:,:) = dcAvg_ded;
    dcAvgStd_ded_group(i,:,:,:) = dcAvgStd_ded;
    
    % Plot prepro HRF results as array map if desired. Make sure you parse the 2D version of the array.
    if plotFlag
        figure;
        DOTHUB_LUMOplotArray(dcAvg_reg,tHRF,prepro.SD2D);
    end
    
    %Convert back to dod for reconstruction
    %dodRecon = DOTHUB_hmrConc2OD(dcAvg_reg,SD3D,[6 6]);
    %tDOD = tHRF;
    
    % Use code snippet from DOTHUB_writePREPRO to define contents of logs:
    %[pathstr, name, ~] = fileparts(nirsFileName);
    %ds = datestr(now,'yyyymmDDHHMMSS');
    %preproFileName = fullfile(pathstr,[name '.prepro']);
    %logData(1,:) = {'Created on: '; ds};
    %logData(2,:) = {'Derived from data: ', nirsFileName};
    %logData(3,:) = {'Pre-processed using:', mfilename('fullpath')};
    %[prepro, preproFileName] = DOTHUB_writePREPRO(preproFileName,logData,dodRecon,tDOD,SD3D,nirs.s,dcAvg,dcAvgStd,tHRF,nirs.CondNames,SD2D);
   
end

%% Do some group averaging here, don't forget to try regressed and non regressed



%%
teststop = 1;

dists = DOTHUB_getSDdists(SD3D);
filt10mm = dists > 7.5 & dists < 12.5;
filt20mm = dists > 17.5 & dists < 22.5;
filt30mm = dists > 27.5 & dists < 32.5;
%% Non-regressed

tmp30 = dcAvg_group(:,:,:,filt30mm);
Av30mm = squeeze(mean(tmp30,4));
newAv30mm = squeeze(mean(Av30mm,1));
tmpStd30 = dcAvgStd_group(:,:,:,filt30mm);
Av30mmStd = squeeze(mean(tmpStd30,4));
newAv30mmStd = squeeze(mean(Av30mmStd,1));
SEM_30mm = newAv30mmStd/sqrt(i);
figure(1);
DOTHUB_plotHRF30(tHRF,newAv30mm,0.5,30)
saveas(gcf,'NonReg_30mm_anag.png');

tmp20 = dcAvg_group(:,:,:,filt20mm);
Av20mm = squeeze(mean(tmp20,4));
newAv20mm = squeeze(mean(Av20mm,1));
tmpStd20 = dcAvgStd_group(:,:,:,filt20mm);
Av20mmStd = squeeze(mean(tmpStd20,4));
newAv20mmStd = squeeze(mean(Av20mmStd,1));
SEM_20mm = newAv20mmStd/sqrt(i);
figure(2);
DOTHUB_plotHRF20(tHRF,newAv20mm,0.5,30)
saveas(gcf,'NonReg_20mm_anag.png');

tmp10 = dcAvg_group(:,:,:,filt10mm);
Av10mm = squeeze(mean(tmp10,4));
newAv10mm = squeeze(mean(Av10mm,1));
tmpStd10 = dcAvgStd_group(:,:,:,filt10mm);
Av10mmStd = squeeze(mean(tmpStd10,4));
newAv10mmStd = squeeze(mean(Av10mmStd,1));
SEM_10mm = newAv10mmStd/sqrt(i);
figure(3);
DOTHUB_plotHRF10(tHRF,newAv10mm,0.5,30)
saveas(gcf,'NonReg_10mm_anag.png');

%% Regressed
tmp30_reg = dcAvg_reg_group(:,:,:,filt30mm);
Av30mm_reg = squeeze(mean(tmp30_reg,4));
newAv30mm_reg = squeeze(mean(Av30mm_reg,1));
tmpStd30_reg = dcAvgStd_reg_group(:,:,:,filt30mm);
Av30mmStd_reg = squeeze(mean(tmpStd30_reg,4));
newAv30mmStd_reg = squeeze(mean(Av30mmStd_reg,1));
SEM_30mm_reg = newAv30mmStd_reg/sqrt(i);
figure(4);
DOTHUB_plotHRF30(tHRF,newAv30mm_reg,0.5,30)
saveas(gcf,'Reg_30mm_anag.png');

tmp20_reg = dcAvg_reg_group(:,:,:,filt20mm);
Av20mm_reg = squeeze(mean(tmp20_reg,4));
newAv20mm_reg = squeeze(mean(Av20mm_reg,1));
tmpStd20_reg = dcAvgStd_reg_group(:,:,:,filt20mm);
Av20mmStd_reg = squeeze(mean(tmpStd20_reg,4));
newAv20mmStd_reg = squeeze(mean(Av20mmStd_reg,1));
SEM_20mm_reg = newAv20mmStd_reg/sqrt(i);
figure(5);
DOTHUB_plotHRF20(tHRF,newAv20mm_reg,0.5,30)
saveas(gcf,'Reg_20mm_anag.png');

tmp10_reg = dcAvg_reg_group(:,:,:,filt10mm);
Av10mm_reg = squeeze(mean(tmp10_reg,4));
newAv10mm_reg = squeeze(mean(Av10mm_reg,1));
tmpStd10_reg = dcAvgStd_reg_group(:,:,:,filt10mm);
Av10mmStd_reg = squeeze(mean(tmpStd10_reg,4));
newAv10mmStd_reg = squeeze(mean(Av10mmStd_reg,1));
SEM_10mm_reg = newAv10mmStd_reg/sqrt(i);
figure(6);
DOTHUB_plotHRF10(tHRF,newAv10mm_reg,0.5,30)
saveas(gcf,'Reg_10mm_anag.png');

%% Regressed - Method 2

tmp30_ded = dcAvg_ded_group(:,:,:,filt30mm);
Av30mm_ded = squeeze(mean(tmp30_ded,4));
newAv30mm_ded = squeeze(mean(Av30mm_ded,1));
tmpStd30_ded = dcAvgStd_ded_group(:,:,:,filt30mm);
Av30mmStd_ded = squeeze(mean(tmpStd30_ded,4));
newAv30mmStd_ded = squeeze(mean(Av30mmStd_ded,1));
SEM_30mm_ded = newAv30mmStd_ded/sqrt(i);
figure(7);
DOTHUB_plotHRF30(tHRF,newAv30mm_ded,SEM_30mm_ded,30)
saveas(gcf,'Ded_30mm_anag.png');

tmp20_ded = dcAvg_ded_group(:,:,:,filt20mm);
Av20mm_ded = squeeze(mean(tmp20_ded,4));
newAv20mm_ded = squeeze(mean(Av20mm_ded,1));
tmpStd20_ded = dcAvgStd_ded_group(:,:,:,filt20mm);
Av20mmStd_ded = squeeze(mean(tmpStd20_ded,4));
newAv20mmStd_ded = squeeze(mean(Av20mmStd_ded,1));
SEM_20mm_ded = newAv20mmStd_ded/sqrt(i);
figure(8);
DOTHUB_plotHRF20(tHRF,newAv20mm_ded,SEM_20mm_ded,30)
saveas(gcf,'Ded_20mm_anag.png');

tmp10_ded = dcAvg_ded_group(:,:,:,filt10mm);
Av10mm_ded = squeeze(mean(tmp10_ded,4));
newAv10mm_ded = squeeze(mean(Av10mm_ded,1));
tmpStd10_ded = dcAvgStd_ded_group(:,:,:,filt10mm);
Av10mmStd_ded = squeeze(mean(tmpStd10_ded,4));
newAv10mmStd_ded = squeeze(mean(Av10mmStd_ded,1));
SEM_10mm_ded = newAv10mmStd_ded/sqrt(i);
figure(9);
DOTHUB_plotHRF10(tHRF,newAv10mm_ded,SEM_10mm_ded,30)
saveas(gcf,'Ded_10mm_anag.png');

%%

% one_participant_plot1 = Av10mm(5,:,:);
% one_par_1 = squeeze(one_participant_plot1);
% figure(4);
% DOTHUB_plotHRF30_onepar(tHRF,one_par_1,SEM_10mm,30)
%saveas(gcf,'NonReg_30mm_par1_stroop.png');

% one_participant_plot2 = Av30mm(2,:,:);
% one_par_2 = squeeze(one_participant_plot2);
% figure(5);
% DOTHUB_plotHRF30_twopar(tHRF,one_par_2,SEM_30mm,30)
% %saveas(gcf,'NonReg_30mm_par2_stroop.png');
% 
% one_participant_plot3 = Av30mm(3,:,:);
% one_par_3 = squeeze(one_participant_plot3);
% figure(6);
% DOTHUB_plotHRF30_onepar(tHRF,one_par_3,SEM_30mm,30)
% %saveas(gcf,'NonReg_30mm_par1_stroop.png');
% 
% one_participant_plot4 = Av30mm(4,:,:);
% one_par_4 = squeeze(one_participant_plot4);
% figure(7);
% DOTHUB_plotHRF30_twopar(tHRF,one_par_4,SEM_30mm,30)
% %saveas(gcf,'NonReg_30mm_par2_stroop.png');
