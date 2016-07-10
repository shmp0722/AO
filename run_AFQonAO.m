function run_AFQonAO

%% Run AFQ_run for these subjects
% set directories

[AFQdata, subs] = SubJect;

%
AMD_Ctl = 9:20;
JMD_Ctl = 38:49;
Ctl     = [24:26,28];
AO      = 22;

SubNumber = [AO,Ctl,AMD_Ctl,JMD_Ctl];
%% Make directory structure for each subject
for ii = 1:length(SubNumber)
    sub_dirs{ii} = fullfile(AFQdata, subs{SubNumber(ii)},'dwi_1st');
end

list = sub_dirs'


% Subject grouping is a little bit funny because afq only takes two groups
% but we have 3. For now we will divide it up this way but we can do more
% later
a = zeros(1,length(SubNumber));
a(1,1) = 1; 
sub_group = a;
% sub_group = [1,0];

% using spm8
addpath(genpath('/home/ganka/bin/spm8'))

% Now create and afq structure
afq = AFQ_Create('sub_dirs', sub_dirs, 'sub_group', sub_group, 'clip2rois',1);
% if you would like to use ants for normalization
% afq = AFQ_Create('sub_dirs', sub_dirs, 'sub_group', sub_group, 'clip2rois', 0,'normalization','ants');

% % To have afq overwrite the old fibers
% afq = AFQ_set(afq,'overwritesegmentation');
% afq = AFQ_set(afq,'overwritecleaning');

% % afq.params.cutoff=[5 95];
% afq.params.outdir = ...
%     fullfile(AFQdata,'/AFQ_results/6LHON_9JMD_8Ctl');
afq.params.outdir='/media/HDPC-UT/dMRI_data/Results/AO';
afq.params.outname = 'afq_29subs.mat';
% afq.params.run_mode = 'mrtrix';

%% test run
afq = AFQ_Create('run_mode','test', 'sub_dirs', sub_dirs, 'sub_group', sub_group); 

%% Run AFQ on these subjects
afq = AFQ_runAO(sub_dirs, sub_group, afq);

%%
afq.params.computenorms =0;

afq = AFQ_SegmentCallosum(afq);
return

%% Next step to see the results
% run 




%% add OR to afq atructure

fName = {'L-OT','R-OT','L-OR','R-OR','LOR0-3','ROR0-3','LOR15-30','ROR15-30'...
    'LOR30-90','ROR30-90'};

roiName1 = {'85_Optic-Chiasm','85_Optic-Chiasm','Lt-LGN4','Rt-LGN4',...
    'Lt-LGN4','Rt-LGN4','Lt-LGN4','Rt-LGN4','Lt-LGN4','Rt-LGN4'};
roiName2 = {'Lt-LGN4','Rt-LGN4','lh_V1_smooth3mm_Half','rh_V1_smooth3mm_Half',...
    'lh_Ecc0to3','rh_Ecc0to3','lh_Ecc15to30','rh_Ecc15to30',...
    'lh_Ecc30to90','rh_Ecc30to90'}

for mm = 1:length(fName)
% L_OT
fgName =fName{mm};
roi1Name =roiName1{mm};
roi2Name = roiName2{mm};
cleanFiber = 0;
computeVals = 1;
showFibers = 0;
segname = [];
overwrite = 1;

afq = AFQ_AddNewFiberGroup(afq, fgName, roi1Name, roi2Name, cleanFiber, computeVals,...
    showFibers,segname,overwrite);

%% R_VOF
fgName = 'R_VOF.pdb';
roi1Name = afq.roi1names{1,20};
roi2Name = afq.roi2names{1,20};
afq = AFQ_AddNewFiberGroup(afq, fgName, roi1Name, roi2Name, 0, 1);




