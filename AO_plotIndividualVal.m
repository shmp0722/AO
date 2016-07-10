function AO_plotIndividualVal

%
%
% make cool individual plot
%
%
%
%
%


%% load afq

% load /sni-storage/wandell/biac2/wandell/data/DWI-Tamagawa-Japan2/RP/afq_Whole_8RP_25Normal_02202015_OTOR.mat
load /media/HDPC-UT/dMRI_data/RP/afq_Whole_8RP_25Normal_02202015_OTOR.mat


%% load TP 

load /media/HDPC-UT/dMRI_data/AO-1-YS-20150531/dwi_1st/fibers/TractProfile_cleanD5_L4.mat


%% Which nodes to analyze
if notDefined('nodes')
    nodes = 21:80;
end

% Define output directory
    outdir = '/sni-storage/wandell/biac2/wandell/data/DWI-Tamagawa-Japan2/RP';

if ~iscell(outdir)
    out_dir{1} = outdir; outdir = out_dir;
end

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
nfg = AFQ_get(afq,'nfg');
% nfg = 20;

fgNames = AFQ_get(afq,'fgnames');

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
%     vout = fullfile(outdir{s},valname{v});
%     if ~exist(vout,'dir')
%         mkdir(vout);
%     end
    
    pVals = AFQ_get(afq,'patient data');
    cVals = AFQ_get(afq,'control data');
    
    % Loop over each fiber group
    for ii = 1:nfg
        % Get the values for the patient and compute the mean
        vals_p = pVals(ii).(upper(valname{v}));
        
        % Remove nodes that are not going to be analyzed and only
        % compute for subject #s
        vals_p = vals_p(:,nodes);
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
        
        for jj = 1:sum(afq.sub_group)
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
   
end

return


%% Render 3d Montage
 
fgNames   = AFQ_get(afq,'fgnames');
tube      = 0;
fgColors  = [];
subNums   = [];
plotDims  = [];
numfibers = 50;

SO_AFQ_MakeFiberGroupMontage(afq, fgNames, tube , fgColors, subNums, plotDims, numfibers);

%% Render 3d figure


% Load the fiber group for the patient
s = 1; % subject

fg_p = AFQ_get(afq,'clean fg',s);
% Load up the b0 image for the patient
dt_p = dtiLoadDt6(AFQ_get(afq,'dt6path',s));
b0_p = readFileNifti(dt_p.files.b0);

% AFQ_RotatingFgGif(fg, colors, outfile, im, slice)
% Make an animated gif of a rotating fiber group

    colors = jet(length(fg_p)+4);

% First we render the first fiber tract in a new figure window
lightH = AFQ_RenderFibers(fg_p(1),'color',colors(1,:),'numfibers',50,'newfig',1);

for ii = 2:length(fg_p)
    % Next add the other fiber tracts to the same figure window
    AFQ_RenderFibers(fg_p(ii),'color',colors(ii,:),'numfibers',50,'newfig',0);
end

% view([0 90])


% lets add OT and OR to the ima
fg{1} = fgRead(afq.files.fibers.LOTD4L4_1206{s});
fg{2} = fgRead(afq.files.fibers.ROTD4L4_1206{s});
fg{3} = fgRead(afq.files.fibers.LOR1206_D4L4{s});
fg{4} = fgRead(afq.files.fibers.ROR1206_D4L4{s});

for ii = 1:4
    AFQ_RenderFibers(fg{ii},'color',colors(length(fg_p)+ii,:),'numfibers',200,'newfig',0);
end

% Add an image if one was provided
AFQ_AddImageTo3dPlot(b0_p,[1 0 0]);
AFQ_AddImageTo3dPlot(b0_p,[0 0 -20]);


view([-48  18])

% Delete the light object and put a new light to the right of the camera
delete(lightH);
lightH=camlight('right');
% Turn of the axes
axis('off');
axis('image');
axis('vis3d');

    