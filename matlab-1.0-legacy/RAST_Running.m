%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Code for running analysis on datasets entered using RAST_Loading.m.
% 
% RAST: Rapid Analysis of Slip Traces, original version written by Vincent Gagneur (2022-2025) v1.0
% The release of the code accompanies a paper which describes its main
% functionalities, and how to interpet the data output in more details, as applied to indentation generated slip traces for bcc/B2 slip activities,  see V. Gagneur et. al. "automated indentation-based slip trace
% analysis for bcc and B2 plasticity" (expected publication in 2025)
% 
% This version of the code is designed analyse SEM micrographs of BCC
% and/or B2 (ordered-bcc) samples, considering {100} {110} {112} and {123} slip plane types simultaneously (48 planes),
% or fcc samples where only {111} slip is expected (4 planes)
% 
% This code requires the "image processing toolbox" from mathworks (install
% via "apps", "get more apps", and find it there). It is possible that it
% depends on other functions from other apps (which I forgot about), just read the error messages when running the
% code and it should be indicated which matlab toolbox you need to install.
% 
% In section (1), User should define here what analysis parameters they
% use. NB the code will always output stats for varying valuest of ST (see
% paper) between 0 and 20.
%
% User can use the loop (2) to run a range of parameters to test (by default,
% varying values of segmentation size N and identification range IR)
%
% A "maintable" containing the results from all different analysis
% parameters is created, and is stored in the main final folder created during the analysis. NB a
% folder is created for each analysis parameter where pictures and individual
% analysis result are saved.
%
% The excel spreadsheet data is layed out to be easily analysed using PivotChart
% in Excel.
%

close all
clear all

prompt = {'Chose a name for your file. It will create a folder where pictures and results are stored'};
dlgtitle = 'Saving the table';
dims = [1 100];

filename = inputdlg(prompt,dlgtitle,dims); %user defines output folder

matfile="demo_B2_tialmo.mat"; % load previously entered eulers + picture combinations e.g. "demo_B2_tialmo.mat". The .mat is created when using the "RAST_Loading.m" matlab code. Make sure to change the name of the file here to match yours, also should be in the same folder as this script.
load(matfile);

%%%%%%%%%%%%%%%%% (1) analysis parameter setting (for pixsize N and error IR, see next loop)

degrot=0; % Default 0 degrees, Anticlockwise rotation to apply to detected line orientations, use if you have a known sample rotation between SEM and EBSD imaging (e.g. if you changed microscope)

nsmooth=10; % Default 10 smoothings, How many smoothing are reapplied to the pictures prior to line detection, if very dirty sample or poorly contrasted microstructure gets detected along slip traces, you may want to increase this to e.g. 20

nei_deg=3.6; % Default 3.6 degrees, Angular difference between neighbouring line orientations to consider them as "similarly oriented" for scoring.

numlines=3; % Default 3 lines, amount of lines detected in each segments.

scorefigs=[0,10]; % default [0,5,10,15], list of score thresholds for which to plot figures (NB, stats for each ST between 0 to 20 will be in the output excel file regardless of this value).

crystal=1; % Crystal geometry, 1 = bcc, 2 = fcc 

%pixsize=100; %Default 100 px, Segmentation size N. Define in loop below by default.

%error=10; %Default 10 degrees, Identification Range IR. Defined in loop below by default.

%%%%%%%%%%%%%%%

maintable=table(); % will store all data
mainfolder=filename{1};
mkdir(filename{1}); %create mainfolder

%%%%% (2) Looping through different code parameters (by default varying N and IR)

for pixsize=100:20:140 % default 100:20:100 segmentation sizes to test e.g. 100:20:140 will test N=100, 120 and 140 px for each identification range IR. For simple analysis, recommended pixsize>90 e.g. = 100

    for error=10:1:10 % default 10:1:10 identification range to test e.g., 1:1:20 will test IR = 1,2,3,...,20 while 6:2:10 will test IR = 6, 8 and 10 for each segmentation size N. For simple analysis use IR = 10

    foldername=strcat(filename,'_',num2str(pixsize),'_N_',num2str(error),'_IR');
    tempT=RAST_function(matfile,mainfolder,foldername,error,pixsize,degrot,nsmooth,numlines,nei_deg,scorefigs,crystal); % calls the "function" version of the main code.
    maintable=[maintable;tempT]; % autorunthrough also saves the information of which "error" (IR) and pxsize (N) that was entered for easy plotting in pivot charts.

    end
end

%%%%%%%%%%%%%%%%

%%% Now making the final spreadsheet

path=strcat(pwd,'\',filename{1});

imname=sprintf('%s_Tablemaster.xlsx', filename{1});
locname=fullfile(path, imname);

writetable(maintable,locname,'Sheet',1,'Range','A1')

close all