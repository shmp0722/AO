
SavePath = '/media/HDPC-UT/dMRI_data/Results/AO/All controls';

Vals = {'fa','md','ad','rd'};
for kk = 1:4;
    for ii = 1:10;
        plotDiffuisonMeasureAO(Vals{kk},ii,SavePath)
    end
end


%% 
SavePath = '/media/HDPC-UT/dMRI_data/Results/AO/Young control';

Vals = {'fa','md','ad','rd'};
for kk = 1:4;
    for ii = 1:10;
        plotDiffuisonMeasureAO_young(Vals{kk},ii,SavePath)
    end
end
