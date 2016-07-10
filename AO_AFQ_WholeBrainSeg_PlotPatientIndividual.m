function AO_AFQ_WholeBrainSeg_PlotPatientIndividual
% Plot patient data against controls
%
%
% Example:
%
% load /biac4/wandell/data/WH/analysis/AFQ_WestonHavens_Full.mat
% afq_controls = afq;
% load /biac4/wandell/data/WH/kalanit/PS/AFQ_PS.mat
% afq_patient = afq;
% AFQ_PlotPatientMeans(afq_patient,afq_controls,'T1_map_lsq_2DTI',[],'age', [53 73])
% AFQ_PlotPatientMeans(afq_patient,afq_controls,'fa',[],'age', [53 73])
% AFQ_PlotPatientMeans(afq_patient,afq_controls,'md',[],'age', [53 73])


%% load afq

% load /sni-storage/wandell/biac2/wandell/data/DWI-Tamagawa-Japan2/RP/afq_Whole_8RP_25Normal_02202015_OTOR.mat
load /media/HDPC-UT/dMRI_data/RP/afq_Whole_8RP_25Normal_02202015_OTOR.mat

afqC = afq;

% load '/media/HDPC-UT/dMRI_data/Results/AO/afq_19-Mar-2016.mat';
load '/media/HDPC-UT/dMRI_data/Results/AO/afq_29subs.mat';


afqP = afq;

clear afq;

%% Which nodes to analyze
if notDefined('nodes')
    nodes = 21:80;
end

% Define output directory
%     outdir = '/media/HDPC-UT/dMRI_data/Results/AO';

% If no value was defined than find all the values that match for controls
% and the patient
if notDefined('valname')
    valname = {'fa' 'md' 'rd' 'ad'};
elseif ischar(valname)
    tmp = valname;clear valname
    valname{1} = tmp;
end
% Get rid of underscores in the valname for the sake of axis labels
for v = 1:length(valname)
    valtitle{v} = valname{v};
    sp = strfind(valname{v},'_');
    if ~isempty(sp)
        valtitle{v}(sp) = ' ';
    end
end
% Get number of fiber groups and their names
nfg = AFQ_get(afqP,'nfg');
% nfg = 28;

fgNames = AFQ_get(afqP,'fgnames');

% Set the views for each fiber group
fgviews = {'leftsag', 'rightsag', 'leftsag', 'rightsag', ...
    'leftsag', 'rightsag', 'leftsag', 'rightsag', 'axial', 'axial',...
    'leftsag', 'rightsag', 'leftsag', 'rightsag',  'leftsag', 'rightsag'...
    'leftsag', 'rightsag', 'leftsag', 'rightsag'};
% Slices to add to the rendering
slices = [-5 0 0; 5 0 0; -5 0 0; 5 0 0; -5 0 0; 5 0 0; -5 0 0; 5 0 0;...
    0 0 -5; 0 0 -5; -5 0 0; 5 0 0; -5 0 0; 5 0 0; -5 0 0; 5 0 0; -5 0 0; 5 0 0; -5 0 0; 5 0 0];

% Set the colormap and color range for the renderings
cmap = AFQ_colormap('bgr');
% cmap = lines(256);

% cmap =colormap;
crange = [-4 4];

% Make an output directory for this subject if there isn't one
% if ~exist(outdir{s},'dir')
%     mkdir(outdir{s});
% end
% fprintf('\nImages will be saved to %s\n',outdir);

%% Loop over the different values
mrvNewGraphWin;

for v = 1:length(valname)
    % Open a new figure window for the mean plot
    subplot(2,2,v); hold('on');
    %     % Make an output directory if it does not exist
    %     vout = fullfile(outdir,valname{v});
    %     if ~exist(vout,'dir')
    %         mkdir(vout);
    %     end
    
    %     pVals = AFQ_get(afqP,'patient data');
    cVals = AFQ_get(afqC,'control data');
    
    % Loop over each fiber group
    for ii = 1:nfg
        % Get the values for the patient and compute the mean
        %         vals_p = pVals(ii).(upper(valname{v}));
        vals_p = afqP.vals.(valname{v}){ii}(1,:);
        %         pVals(ii).(upper(valname{v}));
        
        % Remove nodes that are not going to be analyzed and only
        % compute for subject #s
        vals_p = vals_p(1,nodes);
        %         vals_pm = nanmean(vals_p(:));
        
        % Get the value for each control and compute the mean
        vals_c = cVals(ii).(upper(valname{v}));
        vals_c = vals_c(:,nodes);
        vals_cm = nanmean(vals_c,2);
        
        % Compute control group mean and sd
        m = nanmean(vals_cm);
        sd = nanstd(vals_cm);
        
        % Plot control group means and sd
        x = [ii-.2 ii+.2 ii+.2 ii-.2 ii-.2];
        y1 = [m-sd m-sd m+sd m+sd m-sd];
        y2 = [m-2*sd m-2*sd m+2*sd m+2*sd m-2*sd];
        fill(x,y2, [.6 .6 .6],'edgecolor',[0 0 0]);
        fill(x,y1,[.4 .4 .4] ,'edgecolor',[0 0 0]);
        
        % plot individual means
        %         c = lines(8);
        
        for jj = 1:sum(afqP.sub_group)
            vals_cur = vals_p(jj,:);
            m_curr   = nanmean(vals_cur);
            % Define the color of the point for the fiber group based on its zscore
            tractcol = vals2colormap((m_curr - m)./sd,cmap,crange);
            
            % Plot patient
            plot(ii, m_curr,'ko', 'markerfacecolor',tractcol,'MarkerSize',6);
        end
    end
    
    % make fgnames shorter
    newfgNames = {'l-TR','r-TR','l-C','r-C','l-CC','r-CC','l-CH','r-CH','CFMa',...
        'CFMi','l-IFOF','r-IFOF','l-ILF','r-ILF','l-SLF','r-SLF','l-U','r-U',...
        'l-A','r-A'};
    
    %     set(gca,'xtick',1:nfg,'xticklabel',newfgNames,'xlim',[0 nfg+1],'fontname','times','fontsize',11);
    set(gca,'xtick',1:nfg,'xticklabel',newfgNames,'xlim',[0 nfg+1],'fontname','times','fontsize',11);
    set(gca, 'XTickLabelRotation',90)
    ylabel(upper(valtitle{v}));
    
    h = colorbar('AxisLocation','out');
    h.Label.String = 'z score';
    
    
    
end

return

%% Ok, let's take just callosal fibers
% Loop over the different values
mrvNewGraphWin;

for v = 1:length(valname)
    % Open a new figure window for the mean plot
    subplot(2,2,v); hold('on');
    %     % Make an output directory if it does not exist
    %     vout = fullfile(outdir,valname{v});
    %     if ~exist(vout,'dir')
    %         mkdir(vout);
    %     end
    
    %     pVals = AFQ_get(afqP,'patient data');
    cVals = AFQ_get(afqC,'control data');
    
    % Loop over each fiber group
    for ii = 21:nfg
        % Get the values for the patient and compute the mean
        %         vals_p = pVals(ii).(upper(valname{v}));
        vals_p = afqP.vals.(valname{v}){ii}(1,:);
        %         pVals(ii).(upper(valname{v}));
        
        % Remove nodes that are not going to be analyzed and only
        % compute for subject #s
        vals_p = vals_p(1,nodes);
        %         vals_pm = nanmean(vals_p(:));
        
        % Get the value for each control and compute the mean
        vals_c = cVals(ii).(upper(valname{v}));
        vals_c = vals_c(:,nodes);
        vals_cm = nanmean(vals_c,2);
        
        % Compute control group mean and sd
        m = nanmean(vals_cm);
        sd = nanstd(vals_cm);
        
        % Plot control group means and sd
        x = [ii-.2 ii+.2 ii+.2 ii-.2 ii-.2];
        y1 = [m-sd m-sd m+sd m+sd m-sd];
        y2 = [m-2*sd m-2*sd m+2*sd m+2*sd m-2*sd];
        fill(x,y2, [.6 .6 .6],'edgecolor',[0 0 0]);
        fill(x,y1,[.4 .4 .4] ,'edgecolor',[0 0 0]);
        
        % plot individual means
        %         c = lines(8);
        
        for jj = 1:sum(afqP.sub_group)
            vals_cur = vals_p(jj,:);
            m_curr   = nanmean(vals_cur);
            % Define the color of the point for the fiber group based on its zscore
            tractcol = vals2colormap((m_curr - m)./sd,cmap,crange);
            
            % Plot patient
            plot(ii, m_curr,'ko', 'markerfacecolor',tractcol,'MarkerSize',6);
        end
    end
    
    % make fgnames shorter
    newfgNames = AFQ_get(afqP,'fgNames');
    newfgNames = newfgNames(21:28)
    for kk = 1:length(newfgNames)
        newfgNames{kk} (strfind(newfgNames{kk},'_'))= ' ';
    end
    
    %     set(gca,'xtick',1:nfg,'xticklabel',newfgNames,'xlim',[0 nfg+1],'fontname','times','fontsize',11);
    set(gca,'xtick',21:28,'xticklabel',newfgNames,'fontname','times','fontsize',11);
    set(gca, 'XTickLabelRotation',90)
    ylabel(upper(valtitle{v}));
    
    h = colorbar('AxisLocation','out');
    h.Label.String = 'z score';
    
    
    
end

return

% %% Render 3d figure Axial view
% % Load the fiber group for the patient
% s = 1; % subject
%
% fg_p = AFQ_get(afqP,'clean fg',s);
% % Load up the b0 image for the patient
% dt_p = dtiLoadDt6(AFQ_get(afqP,'dt6path',s));
% b0_p = readFileNifti(dt_p.files.b0);
%
% % AFQ_RotatingFgGif(fg, colors, outfile, im, slice)
% % Make an animated gif of a rotating fiber group
%
% colors = jet(length(fg_p)+4);
%
% % First we render the first fiber tract in a new figure window
% figure; hold on;
% lightH = AFQ_RenderFibers(fg_p(1),'color',colors(1,:),'numfibers',50,'newfig',0);
%
% for ii = 2:length(fg_p)
%     % Next add the other fiber tracts to the same figure window
%     AFQ_RenderFibers(fg_p(ii),'color',colors(ii,:),'numfibers',50,'newfig',0);
% end
%
% % Add callossal fibers
%
% % Add an image if one was provided
% % AFQ_AddImageTo3dPlot(b0_p,[1 0 0]);
% AFQ_AddImageTo3dPlot(b0_p,[0 0 -20]);
%
% view([0 90])
%
% % Turn of the axes
% axis('off');
% axis('image');
% axis('vis3d');
%
% %% Delete the light object and put a new light to the right of the camera
% delete(lightH);
% lightH=camlight('right');
%
% %% side view
% % First we render the first fiber tract in a new figure window
% figure; hold on;
% lightH = AFQ_RenderFibers(fg_p(1),'color',colors(1,:),'numfibers',50,'newfig',0);
%
% for ii = 2:length(fg_p)
%     % Next add the other fiber tracts to the same figure window
%     AFQ_RenderFibers(fg_p(ii),'color',colors(ii,:),'numfibers',50,'newfig',0);
% end
%
% % Add an image if one was provided
% AFQ_AddImageTo3dPlot(b0_p,[1 0 0]);
% % AFQ_AddImageTo3dPlot(b0_p,[0 0 -20]);
%
% view([-90 0])
%
% % Turn of the axes
% axis('off');
% axis('image');
% axis('vis3d');
%
% lightH=camlight('left');

%% Render AO's callosal fibers
figure; hold on;
% get coallosal fiber names
FGN = fieldnames(afqP.files.fibers);

% fg_p = AFQ_get(afqP,'clean fg',s);

% First we render the first fiber tract in a new figure window
% fg = fgRead(afqP.files.fibers.(FGN{5}){s});
% lightH = AFQ_RenderFibers(fg,'color',colors(1,:),'numfibers',50,'newfig',0);

colors = jet(8);

for  kk = 3:10
    fg = fgRead(afqP.files.fibers.(FGN{kk*2}){1});
    AFQ_RenderFibers(fg,'color',colors(kk-2,:),'numfibers',100,'newfig',0);
end

% % Make a CC roi
% roiDir = '/media/HDPC-UT/dMRI_data/AO-1-YS-20150531/ROIs';
% rois = dir(fullfile(roiDir,'*_CC_*mat'));
%
% roi = dtiNewRoi();
% for jj = 1:5
%     cur_ROI = dtiReadRoi(fullfile('/media/HDPC-UT/dMRI_data/AO-1-YS-20150531/ROIs',rois(jj).name))
%
%     roi = dtiMergeROIs(roi, cur_ROI)
% end
% AFQ_RenderRoi(roi)

% Add an image if one was provided
T1 = afqP.files.dt6{1};
dt6 = dtiLoadDt6(T1);
T1 = niftiRead(dt6.files.t1);

% AFQ_AddImageTo3dPlot(b0_p,[-5 0 0]); % b0 imgage
AFQ_AddImageTo3dPlot(T1,[5 0 0]); % t1 image
% AFQ_AddImageTo3dPlot(T1,[0 0 -20]); % t1 image


% Turn of the axes
axis('off');
axis('image');
axis('vis3d');

view([-90 0])
% view([0 90])
% view([-180 45])

% view([-77 22])
% camlight left

title 'AO'
%% Render a control callosal fibers

for ll =3;% 2:28
    figure; hold on;
    % get coallosal fiber names
    FGN = fieldnames(afqP.files.fibers);
    
    % fg_p = AFQ_get(afqP,'clean fg',s);
    
    % First we render the first fiber tract in a new figure window
    % fg = fgRead(afqP.files.fibers.(FGN{5}){s});
    % lightH = AFQ_RenderFibers(fg,'color',colors(1,:),'numfibers',50,'newfig',0);
    
    colors = jet(8);
    
    for  kk = 3:10
        fg = fgRead(afqP.files.fibers.(FGN{kk*2}){ll}); % sub_dirs{26} = YM
        AFQ_RenderFibers(fg,'color',colors(kk-2,:),'numfibers',100,'newfig',0);
    end
    
    %
    
    T1 = afqP.files.dt6{ll};
    T1 = dtiLoadDt6(T1);
    T1 = niftiRead(T1.files.t1);
    
    % Add an image if one was provided
    % AFQ_AddImageTo3dPlot(b0_p,[-5 0 0]); % b0 imgage
    
    AFQ_AddImageTo3dPlot(T1,[10 0 0]); % t1 image
    %     AFQ_AddImageTo3dPlot(T1,[0 0 -20]);
    % Turn of the axes
    axis('off');
    axis('image');
    axis('vis3d');
    
    % camlight left
    hold off;
    
    title(sprintf('Ctl-%d',ll))
    %     view([0 90])
    view([-90 0])
    
    
end
%% Step 5: Render Tract Profiles

% Render the Tract FA Profile for the left arcuate fasciculus. When the
% argument 'dt' is passed in follwed by a variable containing the dt6
% structure then the tract profile is added to the plot. The colormap
% denotes the FA value at each point along the tract core.
figure;
for ii =1:20
    % subplot(4,5,ii)
    
    AFQ_RenderFibers(fg_p(ii),'dt',dt_p,'numfibers',50,'newfig',1);
    title(fgNames{ii})
    xlabel('FA')
    grid off
    axis off
    
end


%% Run a zscore to compare FA along each fiber tract for patients vs. controls

% Get the value for each control and compute the mean
vals_c = cVals(ii).(upper(valname{v}));
vals_c = vals_c(:,nodes);
vals_cm = nanmean(vals_c,2);

% Compute control group mean and sd
m = nanmean(vals_cm);
sd = nanstd(vals_cm);

vals_cur = vals_p(jj,:);
m_curr   = nanmean(vals_cur);
% Define the color of the point for the fiber group based on its zscore
tractcol = vals2colormap((m_curr - m)./sd,cmap,crange);

% Loop over all 20 fiber groups
for jj = 1:nfg
    % Run an independent samples t-test comparing FA values between the
    % groups at each point along each tract
    [h(jj,:),p(jj,:),~,Tstats(jj)] = zscore(afqP.patient_data(jj).FA,afqC.control_data(jj).FA);
end

%% Render figure showing OT OR
% Load ACH

load /media/HDPC-UT/dMRI_data/Results/ACH_0402.mat

%% load OT, OR

otDir = '/media/HDPC-UT/dMRI_data/AO-1-YS-20150531/dwi_1st/fibers/conTrack/OT_5K';
orDir = '/media/HDPC-UT/dMRI_data/AO-1-YS-20150531/dwi_1st/fibers/conTrack/OR_100K';

% load lot rot
lot = fgRead(fullfile(otDir,[ACH{22,1}.name,'.mat']));

rot =fgRead(fullfile(otDir,...
    'fg_OT_5K_85_Optic-Chiasm_Rt-LGN4_2015-07-13_18.48.12-2_Left-Cerebral-White-Matter_Ctrk100_cleaned_deplicated.mat'));

% load OR
lor = fgRead(fullfile(orDir,[ACH{22,3}.name,'.pdb']));
ror = fgRead(fullfile(orDir,[ACH{22,4}.name,'.pdb']));


% Render AO's OT OR
figure; hold on;
% get coallosal fiber names
% FGN = fieldnames(afqP.files.fibers);

colors = jet(2);

% for  kk = 3:10
AFQ_RenderFibers(lot,'color',colors(1,:),'numfibers',100,'newfig',0);
AFQ_RenderFibers(rot,'color',colors(1,:),'numfibers',100,'newfig',0);

AFQ_RenderFibers(lor,'color',colors(2,:),'numfibers',100,'newfig',0);
AFQ_RenderFibers(ror,'color',colors(2,:),'numfibers',100,'newfig',0);

% Add an image if one was provided
T1 = afqP.files.dt6{1};
dt6 = dtiLoadDt6(T1);
T1 = niftiRead(dt6.files.t1);

% AFQ_AddImageTo3dPlot(b0_p,[-5 0 0]); % b0 imgage
% AFQ_AddImageTo3dPlot(T1,[1 0 0]); % t1 image
AFQ_AddImageTo3dPlot(T1,[0 0 -15]); % t1 image


% Turn of the axes
axis('off');
axis('image');
axis('vis3d');

% view([-90 0])
% view([0 90])
view([-142 32])
% view([-180 45])

% view([-77 22])
% camlight left

title 'AO'

%% Render Ctl 13
ll =3;% ACH{24}
%     figure; hold on;

otDir = '/media/HDPC-UT/dMRI_data/Ctl-13-MW-20140313-dMRI-Anatomy/dwi_1st/fibers/conTrack/OT_5K';
orDir = '/media/HDPC-UT/dMRI_data/Ctl-13-MW-20140313-dMRI-Anatomy/dwi_1st/fibers/conTrack/OR_100K';

% load lot rot
lot = fgRead(fullfile(otDir,[ACH{24,1}.name,'.mat']));
rot = fgRead(fullfile(otDir,[ACH{24,2}.name,'.mat']));


% load OR
lor = fgRead(fullfile(orDir,[ACH{24,3}.name,'.pdb']));
ror = fgRead(fullfile(orDir,[ACH{24,4}.name,'.pdb']));


% Render AO's OT OR
figure; hold on;
% get coallosal fiber names
% FGN = fieldnames(afqP.files.fibers);

colors = jet(2);

% for  kk = 3:10
AFQ_RenderFibers(lot,'color',colors(1,:),'numfibers',100,'newfig',0);
AFQ_RenderFibers(rot,'color',colors(1,:),'numfibers',100,'newfig',0);

AFQ_RenderFibers(lor,'color',colors(2,:),'numfibers',100,'newfig',0);
AFQ_RenderFibers(ror,'color',colors(2,:),'numfibers',100,'newfig',0);

% Get dt6 and t1
T1 = afqP.files.dt6{ll};
dt6 = dtiLoadDt6(T1);
T1 = niftiRead(dt6.files.t1);

% Add an image if one was provided
% AFQ_AddImageTo3dPlot(b0_p,[-5 0 0]); % b0 imgage

% AFQ_AddImageTo3dPlot(T1,[10 0 0]); % t1 image
    AFQ_AddImageTo3dPlot(T1,[0 0 -20]);
% Turn of the axes
axis('off');
axis('image');
axis('vis3d');

% camlight left
hold off;

title(sprintf('Ctl-%d',ll))
    view([0 90])
% view([-90 0])

%% Render AO's callosal fibers
% figure; hold on;
% get coallosal fiber names
FGN = fieldnames(afqP.files.fibers);

% fg_p = AFQ_get(afqP,'clean fg',s);

% First we render the first fiber tract in a new figure window
% fg = fgRead(afqP.files.fibers.(FGN{5}){s});
% lightH = AFQ_RenderFibers(fg,'color',colors(1,:),'numfibers',50,'newfig',0);

colors = jet(8);

% Add an image if one was provided
T1 = afqP.files.dt6{1};
dt6 = dtiLoadDt6(T1);
T1 = niftiRead(dt6.files.t1);


for  kk = 3:10
    fg = fgRead(afqP.files.fibers.(FGN{kk*2}){s});
    AFQ_RenderFibers(fg,'dt',dt6,'color',colors(kk-2,:),'numfibers',100,'newfig',1);
end

% % Make a CC roi
% roiDir = '/media/HDPC-UT/dMRI_data/AO-1-YS-20150531/ROIs';
% rois = dir(fullfile(roiDir,'*_CC_*mat'));
%
% roi = dtiNewRoi();
% for jj = 1:5
%     cur_ROI = dtiReadRoi(fullfile('/media/HDPC-UT/dMRI_data/AO-1-YS-20150531/ROIs',rois(jj).name))
%
%     roi = dtiMergeROIs(roi, cur_ROI)
% end
% AFQ_RenderRoi(roi)




% AFQ_AddImageTo3dPlot(b0_p,[-5 0 0]); % b0 imgage
AFQ_AddImageTo3dPlot(T1,[5 0 0]); % t1 image
% AFQ_AddImageTo3dPlot(T1,[0 0 -20]); % t1 image


% Turn of the axes
axis('off');
axis('image');
axis('vis3d');

view([-90 0])
% view([0 90])
% view([-180 45])

% view([-77 22])
% camlight left

title 'AO'

%%
% Render the Tract FA Profile for the left corticospinal tract
AFQ_RenderFibers(fg_p(3),'dt',dt_p);

% Render the Tract FA Profile for the left IFOF
AFQ_RenderFibers(fg_p(11),'dt',dt_p);

% Render the Tract FA Profile for the left uncinate
AFQ_RenderFibers(fg_p(17),'dt',dt_p);


end

