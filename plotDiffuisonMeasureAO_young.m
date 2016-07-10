function plotDiffuisonMeasureAO_young(vals,fibID,SavePath)
% Plot figure 5 showing individual FA value along the core of OR and optic tract.
%
% Repository dependencies
%    VISTASOFT
%    AFQ
%    LHON2
%
% SO Vista lab, 2014
% 
% vals = 'fa'; % fa ,md ,ad ,rd 
% fibID = 1; 1:10;
% 
% plotDiffuisonMeasureAO(vals,fibID,1); save
% 
% Shumpei Ogawa 2014

%% Identify the directories and subject types in the study
% The full call can be
[homeDir, subDir, ~, AMD_Ctl, ~, Ctl] = SubJect20160128;


% Load ACH data
TPdata = '/media/HDPC-UT/dMRI_data/Results/ACH_0827.mat';
load(TPdata);

%
if notDefined('vals')
    vals = 'fa';
end

if notDefined('fibID')
    fibID = 1;
end

%% Figure
% indivisual FA value along optic tract
% if fibID< 5,
% take values
fbName = {'L-OT','R-OT','L-OR','R-OR','LOR0-3','ROR0-3','LOR15-30','ROR15-30'...
    'LOR30-90','ROR30-90'};
% package to cnotain
nodes =  length(ACH{22,fibID}.vals.fa);
fa = nan(length(subDir), nodes);
md = fa;
ad = fa;
rd = fa;

% unite values 
for subID = 1:length(ACH);
    if isempty(ACH{subID,fibID});
        fa(subID,:) =nan(1,nodes);
    else
        fa(subID,:) =  ACH{subID,fibID}.vals.fa;
    end;
    
    if isempty(ACH{subID,fibID});
        md(subID,:) =nan(1,nodes);
    else
        md(subID,:) = ACH{subID,fibID}.vals.md;
    end;
    
    if isempty(ACH{subID,fibID});
        rd(subID,:) =nan(1,nodes);
    else
        rd(subID,:) = ACH{subID,fibID}.vals.rd;
    end;
    
    if isempty(ACH{subID,fibID});
        ad(subID,:) =nan(1,nodes);
    else
        ad(subID,:) = ACH{subID,fibID}.vals.ad;
    end;
end

vals = lower(vals);
switch vals
    case 'fa'
        val_C  = fa(Ctl,:);
        val_AC = fa(AMD_Ctl,:);
        AO     = fa(22,:);
%         val_RP = fa(RP,:);
    case 'md'
        val_C  = md(Ctl,:);
        val_AC = md(AMD_Ctl,:);
        AO     = md(22,:);
%         val_RP = md(RP,:);
    case 'ad'
        val_C  = ad(Ctl,:);
        val_AC = ad(AMD_Ctl,:);
        AO     = ad(22,:);

%         val_RP = ad(RP,:);
    case 'rd'
        val_C  = rd(Ctl,:);
        val_AC = rd(AMD_Ctl,:);
        AO     = rd(22,:);

%         val_RP = rd(RP,:);
end

%
CTL_data = [val_C];% [val_C;val_AC];
% RP_data  = val_RP;

% %% ANOVA
% group =2;
% for jj= 1: nodes
%     M = length(Ctl)+length(AMD_Ctl);
%     pac = nan(M,group);
%     pac(:,1)= [val_C(:,jj);val_AC(:,jj)];
%     pac(1:8,2)= val_RP(:,jj);
%     
%     [p(jj),~,stats(jj)] = anova1(pac,[],'off');
%     co = multcompare(stats(jj),'display','off');
%     C{jj}=co;
% end
% Portion =  p<0.05; % where is most effected
% 
% Portion = Portion+0;
%% Render individuals and Normal range
h1= figure; hold on;
X = 1:nodes;
c = lines(100);

% % put bars based on ANOVA (p<0.01)
% B= bar(1:nodes,Portion*2.5,1);
% set(B,'EdgeColor','none')

% Control
st = nanstd(CTL_data,1);
m   = nanmean(CTL_data,1);

% render control subjects range
A3 = area(m+2*st);
A1 = area(m+st);
A2 = area(m-st);
A4 = area(m-2*st);

% set color and style
set(A1,'FaceColor',[0.6 0.6 0.6],'linestyle','none')
set(A2,'FaceColor',[0.8 0.8 0.8],'linestyle','none')
set(A3,'FaceColor',[0.8 0.8 0.8],'linestyle','none')
set(A4,'FaceColor',[1 1 1],'linestyle','none')

plot(m,'color',[0 0 0], 'linewidth',3 )

% add individual FA plot
% for k = 22 %1:length(subDir)
    plot(X,AO,'Color',c(3,:),...
        'linewidth',1);
% end
% mRP   = nanmean(RP_data,1);
% plot(X,mRP,'Color',c(3,:) ,'linewidth',3)

T= title(sprintf('%s', fbName{fibID}));
ylabel(upper(vals))
xlabel('Location')

% set axes
a =[10,40];
set(gca,'xlim',a,'xTick',a)
set(gca,'ytickMode','auto')
b =get(gca,'ylim');

% set Y limit 
X = (m+2*st);
A =  round(max(X(10:40))+0.2,1);

Y = (m-2*st);
B =  round(min(Y(10:40))-0.2,1);

b = [B,A];

if b(1)<0;b(1)=0;end;

set(gca,'ylim',b,'yTick',b);

hold off;

% save the figure
if ~isempty(SavePath)
%     saveas(h1,fullfile(SavePath, [vals,'_',T.String]),'pdf')
%     saveas(h1,fullfile(SavePath, [vals,'_',T.String,'.eps']),'psc2')
    saveas(h1,fullfile(SavePath, [vals,'_',T.String]),'bmp')
end

return

%% Optic Tract compare to Ctl
figure; hold on;

% Control
st = nanstd(val_C);
m   = nanmean(val_C,1);

% render control subjects range
A3 = area(m+2*st);
A1 = area(m+st);
A2 = area(m-st);
A4 = area(m-2*st);

% set color and style
set(A1,'FaceColor',[0.6 0.6 0.6],'linestyle','none')
set(A2,'FaceColor',[0.8 0.8 0.8],'linestyle','none')
set(A3,'FaceColor',[0.8 0.8 0.8],'linestyle','none')
set(A4,'FaceColor',[1 1 1],'linestyle','none')

plot(m,'color',[0 0 0], 'linewidth',3 )

% add individual FA plot
for k = 1:length(RP) %1:length(subDir)
    plot(X,RP_data(k,:),'Color',c(3,:),...
        'linewidth',1);
end
m   = nanmean(RP_data,1);
plot(X,m,'Color',c(3,:) ,'linewidth',3)

T = title(sprintf('%s comparing to Ctl', fbName{fibID}));
ylabel(upper(vals))
xlabel('Location')

% set axes
a =[10,40];
set(gca,'xlim',a,'xTick',a)
set(gca,'ytickMode','auto')
b =get(gca,'ylim');
if b(1)<0;b(1)=0;end;

set(gca,'ylim',b,'yTick',b);

hold off;

% Save current figure
if ~isempty(SavePath)
    saveas(h1,fullfile(SavePath, [vals,'_',T.String]),'pdf')
    saveas(h1,fullfile(SavePath, [vals,'_',T.String,'.eps']),'psc2')
    saveas(h1,fullfile(SavePath, [vals,'_',T.String]),'bmp')
    
end

%% Optic Tract compare to AMD_Ctl
figure; hold on;

% Control
st = nanstd(val_AC);
m   = nanmean(val_AC,1);

% render control subjects range
A3 = area(m+2*st);
A1 = area(m+st);
A2 = area(m-st);
A4 = area(m-2*st);

% set color and style
set(A1,'FaceColor',[0.6 0.6 0.6],'linestyle','none')
set(A2,'FaceColor',[0.8 0.8 0.8],'linestyle','none')
set(A3,'FaceColor',[0.8 0.8 0.8],'linestyle','none')
set(A4,'FaceColor',[1 1 1],'linestyle','none')

plot(m,'color',[0 0 0], 'linewidth',3 )

% add individual FA plot
for k = 1:length(RP) %1:length(subDir)
    plot(X,RP_data(k,:),'Color',c(3,:),...
        'linewidth',1);
end
m   = nanmean(RP_data,1);
plot(X,m,'Color',c(3,:) ,'linewidth',3)

T = title(sprintf('%s comparing to AMD_C', fbName{fibID}));
ylabel(upper(vals))
xlabel('Location')

% set axes
a =[10,40];
set(gca,'xlim',a,'xTick',a)
set(gca,'ytickMode','auto')
b =get(gca,'ylim');
if b(1)<0;b(1)=0;end;

set(gca,'ylim',b,'yTick',b);

hold off;


% Save current figure
if ~isempty(SavePath)
    saveas(h1,fullfile(SavePath,[vals,'_',T.String]),'pdf')
    saveas(h1,fullfile(SavePath, [vals,'_',T.String,'.eps']),'psc2')
    saveas(h1,fullfile(SavePath, [vals,'_',T.String]),'bmp')
    
end
return
