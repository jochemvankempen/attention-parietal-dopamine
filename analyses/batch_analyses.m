function batch_analyses()
% 

subjects = {'J','W'};% dopamine data
subjects = {'W'};% dopamine data
% subjects = {'S'}; % control data

%% paths

path_data = ['C:\Jochem\Gratc_PPC_DA_new\data\processed'];
path_data = ['/Users/jochemvankempen/NCL/gratc_DA/processed'];

%% addpaths

addpath(genpath(''))

%% get list of recordings

recordinglist = [];

for isub = 1:length(subjects)
    recs = dir(fullfile(path_data, subjects{isub}, '2*'));
    for irec = 1:length(recs)
        recinfo = load(fullfile(recs(irec).folder, recs(irec).name, 'recinfo.mat'));
        recordinglist = [recordinglist; recinfo.recinfo];
    end
end

%% define jobs

Job.spike_rate_ANOVA = 1;
Job.spike_rate_ANOVA = 1;

%% Run jobs

for irec = 1:height(recordinglist)
    
    recinfo = recordinglist(irec,:);
    
    fprintf('Running analyses for: Subject %s, recording %s\n', recinfo.Subject, recinfo.Date)
    fprintf('--------------------------------------------------------\n')
    
    perform_analyses(recinfo, Job, path_data)
        
end

%% population analysis
batch_gratc_populationAnalysis(subjects)


%% Extract single paradigm from full recording

for isubject = 1:length(subjects)
    
    assert(strcmpi(subjects{isubject},'Stitch'), 'Only run data extraction for Stitch')
    
    % get info
    [list_pen, list_unit] = setList_gratc_DA(subjects{isubject});
    
    
    %%  little hard-coded hack to extract control data from raw files from Stitch
    if 0
        
        NEVpath = fullfile(list_pen.basePath, 'NLX_offline', list_pen.replayPath, 'Events.nev');
        
        for irec = 1:height(list_pen)
            NEV = NLX_LoadNEV(NEVpath{irec},'Full',1,[]);
            eventString = {[list_pen.grcjdru1{irec} ' on'], [list_pen.grcjdru1{irec} ' off']};
            
            [Index,TimeStamps,NEV] = NLX_findEvents(NEV, 'EVENTSTRING', eventString);
            
            assert(length(find(Index))==2, 'Did not find correct number of events')
            
            for ichan = 1:3
                % the same signal extracted using different spike thresholds
                
                % load-filename
                NSEpath = fullfile(list_pen.basePath{irec}, 'NLX_offline', list_pen.replayPath{irec}, sprintf('SE%d.NSE', ichan));
                
                % extract data
                [NSE] = NLX_LoadNSE(NSEpath,'FULL',4,TimeStamps);
                
                % save-filename
                NSEpath_save = fullfile(list_pen.basePath{irec}, 'NLX_offline', list_pen.replayPath{irec}, sprintf('SE%d_%s.NSE', ichan, list_pen.grcjdru1{irec}));
                
                NSE.Path = NSEpath_save;
                NLX_SaveNSE(NSE,0,0);
            end
        end
    end
    
    
    %% Pen analysis
    
    for iparadigm = 1:height(list_pen)
        batch_msacc_penAnalysis(subjects{isubject}, list_pen(iparadigm,:), 1)
    end
    %% unit analysis
    
    for iunit = 5:height(list_unit)
        filename = regexp( list_unit.Filename{iunit}, '(?<=_)[^]+(?=.NSE$)','once','match');
        penIdx = strcmpi(list_pen.grcjdru1, filename);

        batch_msacc_unitAnalysis(subjects{isubject}, list_pen(penIdx,:), list_unit(iunit,:), 1)
    end
    %% population analysis
    
    event = 'NLX_SACCADE_START';
    
    all_significance = zeros(height(list_unit), 3);
    all_rates = zeros(height(list_unit), 2, 4);
    
    for iunit = 1:height(list_unit)
        
        filename = regexp( list_unit.Filename{iunit}, '(?<=_)[^]+(?=.NSE$)','once','match');
        penIdx = strcmpi(list_pen.grcjdru1, filename);

        [readpaths, targetpaths] = setPaths_gratc_DA('Stitch', list_pen(penIdx,:));

        loadfilename = fullfile(targetpaths.pen,  sprintf('PEN%d_%s_sps_cell%d.mat', list_unit.PenID{iunit}, list_pen.grcjdru1{penIdx}, list_unit.CellID{iunit}));
        
        unit_data = load(loadfilename);
        
        % find the relevant event
        event_idx = strcmpi(unit_data.event2align, event);
        
        % collect significance levels
        all_significance(iunit,:) = unit_data.Significance;
    
        % collect spike rates
        all_rates(iunit,:,:) = unit_data.spsRate_mean(:,:,event_idx);
        
        
    end
    
    drug_off_cond = squeeze(all_rates(:,1,:));
    drug_on_cond = squeeze(all_rates(:,2,:));
    
    drug_off = squeeze(mean(all_rates(:,1,:),3));
    drug_on = squeeze(mean(all_rates(:,2,:),3));

    
    % plot
    [fH, fSet] = plotj_initFig('figNum', 3, 'journal', 'Science', 'width', 16, 'height', 12);
    
    ncol=4;
    nrow=2;
    
    % across units
    iplot=5;
    subtightplot(nrow,ncol,[iplot iplot+1],fSet.subplotGap,fSet.subplotMargin,fSet.subplotMargin)
    plotj_initAx(fSet, 'axlabel', 2);
    hold on
    plotj_scatter([drug_off(:), drug_on(:)], 'unityLine', 1, 'markersize', fSet.MarkerSize);
    xlabel({'Control','Firing rate (Hz)'}, 'FontSize', fSet.Fontsize_text)
    ylabel({'Saline','Firing rate (Hz)'}, 'FontSize', fSet.Fontsize_text)
    axis square
    [P, STATS, textOut, fullTextOut] = statj_compareMeans(drug_off(:), drug_on(:), 1, 'rank');
    
    text(3,15,textOut,'FontSize',fSet.Fontsize_text)

    % across units
    iplot=7;
    subtightplot(nrow,ncol,[iplot iplot+1],fSet.subplotGap,fSet.subplotMargin,fSet.subplotMargin)
    plotj_initAx(fSet, 'axlabel', 3);
    hold on
    plotj_scatter([drug_off_cond(:), drug_on_cond(:)], 'unityLine', 1, 'markersize', fSet.MarkerSize);
    xlabel({'Control','Firing rate (Hz)'}, 'FontSize', fSet.Fontsize_text)
    ylabel({'Saline','Firing rate (Hz)'}, 'FontSize', fSet.Fontsize_text)
    axis square
    [P, STATS, textOut, fullTextOut] = statj_compareMeans(drug_off_cond(:), drug_on_cond(:), 1, 'rank');
    
    text(3,15,textOut,'FontSize',fSet.Fontsize_text)

    targetpath = 'C:\Jochem\Gratc_PPC_DA_new\figures\';
    savefigname = sprintf('population_sacc');
    figSave(savefigname, targetpath, {'png', 'svg'});

end

