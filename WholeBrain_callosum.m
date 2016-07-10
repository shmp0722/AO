% WholeBrain_callosum

%%

[homeDir, sub_dirs] = SubJect;



%%
parfor ii =1: length(sub_dirs)
    if ~exist(fullfile(sub_dirs{ii},'fibers/WholeBrainFG.mat'),'file')
        AFQ_TrackAndSegmentOneSub(fullfile(homeDir,sub_dirs{ii},'dwi_1st'))
    end
end

%
parfor ii =1: length(sub_dirs)
    if ~exist(fullfile(sub_dirs{ii},'fibers/callosumFG.mat'),'file')
        AFQ_TrackAndSegmentCallosumOneSub_SO(fullfile(homeDir,sub_dirs{ii},'dwi_1st'))
    end
end
