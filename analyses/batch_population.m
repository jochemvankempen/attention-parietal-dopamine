function batch_population()

clear all
addpath(genpath('./attention-parietal-dopamine'))
addpath(genpath('./plotj'))

%% param
subjects = {'J','W'};
% path_data = ['C:\Jochem\Gratc_PPC_DA_new\data\processed'];
path_data = ['/Users/jochemvankempen/NCL/gratc_DA/processed'];
path_population = regexprep(path_data,'processed','population');

if ~isfolder(path_population)
    mkdir(path_population)
end

%% plot param
label_drug = {'Dopamine','SCH23390'};
label_drug_ext = {'',' (D1-antagonist)'};
label_attention = {'RF','away'};
label_drug_onoff = {'Off','On'};
label_celltype = {'Narrow','Broad'};

markersize = 20; % scatterplot markersize

plot_conventions = 'Science'; % follow figure conventions from this journal

%% get recordinglist
recordinglist = get_recordinglist(subjects, path_data);

%% Get population data

[rate_summary, unitlist] = get_population_data(recordinglist, 'spike_rate_summary', path_data, 'dim');
rate_PSTH = get_population_data(recordinglist, 'spike_rate_PSTH', path_data, 'dim');
rate_ANOVA = get_population_data(recordinglist, 'spike_rate_ANOVA', path_data, 'dim');
rate_ROC = get_population_data(recordinglist, 'spike_rate_ROC', path_data, 'dim');
waveform = get_population_data(recordinglist, 'waveform', path_data);
% pupil_windows = get_population_data(recordinglist, 'pupil_windows', path_data, 'stim');
% pupil_timeseries = get_population_data(recordinglist, 'pupil_timeseries', path_data, 'stim');

RT = get_population_data(recordinglist, 'RT_drug_modulation', path_data, 'stim');

%% Extract param

idx_attention = [1 2 2 1 2 2]; % cond_num mapping to attend RF/away in trialdata
idx_attention_roc = [1 1 2 1 1 2]; % see condition comparisons in spike_rate_ROC.m
num_units = size(waveform.waveform,1);
num_recs = size(RT.RT,1);
% selectivity = get_unit_selectivity(rate_ANOVA.p_anova, 'att+dru')

%% classify cells as narrow or broad spiking

waveform_cutoff = 250;

unit_class = ones(num_units,1);
unit_class(waveform.peak_to_trough_time > waveform_cutoff) = 2;


%% figure 1 - plot cell proportions

% att*drug
interaction_selective = get_unit_selectivity(unitlist, rate_ANOVA.selectivity, 'att*dru');

% visual&att&drug
ABC = get_unit_selectivity(unitlist, rate_ANOVA.selectivity, 'visual&att&dru');

% visually selective
A = get_unit_selectivity(unitlist, rate_ANOVA.selectivity, 'visual');

% attention selective
B = get_unit_selectivity(unitlist, rate_ANOVA.selectivity, 'att');

% visual and attention
AB = get_unit_selectivity(unitlist, rate_ANOVA.selectivity, 'visual&att');

% drug
C = get_unit_selectivity(unitlist, rate_ANOVA.selectivity, 'dru');

% att&drug
BC = get_unit_selectivity(unitlist, rate_ANOVA.selectivity, 'att&dru');

% visual&drug
AC = get_unit_selectivity(unitlist, rate_ANOVA.selectivity, 'visual&dru');



ABC_venn    = length(find( ABC ));
AB_venn     = length(find( AB   & ~(ABC) ));
BC_venn     = length(find( BC   & ~(ABC) ));
AC_venn     = length(find( AC   & ~(ABC) ));

A_venn      = length(find( A    & ~(ABC | AB | AC) ));
B_venn      = length(find( B    & ~(ABC | AB | BC) ));
C_venn      = length(find( C    & ~(ABC | AC | BC) ));




[fH, fSet] = plotj_initFig('width', 10, 'height', 10);


subtightplot(1, 1, 1, fSet.subplotGap, fSet.subplotMargin, fSet.subplotMargin)
plotj_initAx(fSet, 'axlabel', 1);
hold on
bar([length(find(A)), length(find(B)), length(find(C)), length(find(interaction_selective))]/num_units)
ylabel('Proportion of units', 'FontSize', fSet.Fontsize_text)
xlim([0.2 4.8])

labels = {...
    'Visual response',...
    'Attention',...
    'Drug',...
    'Interaction'};

set(gca,'XTickLabel',labels,'XTickLabelRotation',-25, 'FontSize', fSet.Fontsize_text)
axis square
% subtightplot(1, 2, 2, fSet.subplotGap, fSet.subplotMargin, fSet.subplotMargin)
% plotj_initAx(fSet, 'axlabel', 2);
%
% venn_prop = [...
%     A_venn, AB_venn, B_venn... % A, A&B, B
%     BC_venn C_venn ... % B&C, C
%     AC_venn ... %C&A
%     ABC_venn
%     ];
% vennX(venn_prop, 1/100)

% Ap = [A B C];
% Ip = [AB AC BC ABC];
% venn(Ap, Ip)


% text(0,0.8,sprintf('%d units with visual response (%1.1f%%)\n', selective_visual, selective_visual/num_units*100))
% set(gca,'xColor',[1 1 1],'yColor',[1 1 1])


% ylim(ylim2use)
savefigname = sprintf('cellProportion');
figSave(savefigname, readpaths.fig, {'png', 'svg'})

%% Print statistics

fprintf('UNITS\n')
fprintf('\tRecorded %d units\n', height(unitlist))

fprintf('SUBJECTS\n')
M1 = strcmpi(unitlist.Subject, {'W'});
M2 = strcmpi(unitlist.Subject, {'J'});
fprintf('\tRecorded %d/%d units from monkey 1/2\n', length(find(M1)), length(find(M2)))

fprintf('SELECTIVITY\n')
unit_selectivity = get_unit_selectivity(unitlist, rate_ANOVA.selectivity, 'visual');
selective_visual = length(find(unit_selectivity));
fprintf('\t%d units with visual response (%1.1f%%)\n', selective_visual, selective_visual/num_units*100)
unit_selectivity = get_unit_selectivity(unitlist, rate_ANOVA.selectivity, 'att');
selective_attention = length(find(unit_selectivity));
fprintf('\t%d units with attentional modulation (%1.1f%%)\n', selective_attention, selective_attention/num_units*100)
% unit_selectivity = get_unit_selectivity(unitlist, rate_ANOVA.selectivity, 'visual&att');
% selective_visAtt = length(find(unit_selectivity));
%
% unit_selectivity = selectData(unitlist, 'visualDrug', 'NA', allSelectivity);
% selective_visDrug = length(find(unit_selectivity));

fprintf('DRUG\n')
selective_drug = 0;
selective_interaction = 0;
for idrug = 1:length(label_drug)
    
    num_drug_units = strcmpi(unitlist.Drug, label_drug(idrug)); % total units recorded with this drug
    
    unit_selectivity = get_unit_selectivity(unitlist, rate_ANOVA.selectivity, 'dru',{'drug',label_drug(idrug)});
    
    fprintf('\tRecorded %d/%d (%d) units from monkey with %s 1/2\n', length(find(M1 & num_drug_units)), length(find(M2 & num_drug_units)), length(find(num_drug_units)), label_drug{idrug})
    
    selective_drug = selective_drug+length(find(unit_selectivity));
    fprintf('\t%d units modulated by %s (%1.1f%%)\n', length(find(unit_selectivity)), label_drug{idrug}, length(find(unit_selectivity))/length(find(num_drug_units))*100)
    
    unit_selectivity = get_unit_selectivity(unitlist, rate_ANOVA.selectivity, 'att*dru',{'drug',label_drug(idrug)});
    selective_interaction = selective_interaction+length(find(unit_selectivity));
    fprintf('\t%d units interaction with attention %s (%1.1f%%)\n', length(find(unit_selectivity)), label_drug{idrug}, length(find(unit_selectivity))/length(find(num_drug_units))*100)
end
fprintf('\t%d total units modulated by attention  (%1.1f%%)\n', selective_attention, selective_attention/num_units*100)
fprintf('\t%d total units modulated by drug application  (%1.1f%%)\n', selective_drug, selective_drug/num_units*100)
fprintf('\t%d total units interaction attention by drug application  (%1.1f%%)\n', selective_interaction, selective_interaction/num_units*100)

%% figure 2 - plot population histogram and example unit

colors = get_colors('att_drug');

events = {'stim', 'cue', 'dim'};
label_event = {'Stimulus onset', 'Cue onset', 'First-dimming'};
xlabel_event = {{'Time from','stimulus onset (ms)'}, {'Time from','cue onset (ms)'}, {'Time from','first-dimming (ms)'}};
xlim2use = {[-200 750], [-100 750], [-1000 200]};
ylim2use = [25 90];
ncol = length(label_event) + 2;
nrow = length(label_drug);

idx_subplot = {1,2,3,6,7,8,[4 5 9 10]};

[fH, fSet] = plotj_initFig('width', 18, 'height', 11, 'Journal',plot_conventions);
fSet.subplotGap2 = fSet.subplotGap .* [1.3 0.30];

iplot = 0;
for irow = 1:nrow
    
    idx_unit = get_unit_selectivity(unitlist, rate_ANOVA.selectivity, 'att&dru', {'drug',label_drug(irow)});
    
    clear h
    for ievent = 1:length(label_event)
        iplot = iplot+1;
        
        rate_PSTH = get_population_data(recordinglist, 'spike_rate_PSTH', path_data, events{ievent});
        
        times = rate_PSTH.time;
        samples = rate_PSTH.samples;
        
        subtightplot(nrow, ncol, idx_subplot{iplot}, fSet.subplotGap2, fSet.subplotMargin, fSet.subplotMargin)
        if ievent==1 && irow==1
            plotj_initAx(fSet, 'axlabel', irow, 'axlabelDisplacement', [0.07, 0.00]);
        else
            plotj_initAx(fSet);
        end
        
        if ievent==(ceil(length(label_event)/2))
            ht = title([label_drug{irow} label_drug_ext{irow}], 'FontSize', fSet.Fontsize_title);
            %             ht.Position = ht.Position + [0.05 0 0];
        end
        if (ievent~=1)
            axh = get(gca);
            axh.YTick = [];
            axh.YAxis.Visible = 'off'; % remove y-axis
            %             axes('Color','none','YColor','none')
            %             set(gca,'xtick',[])
        end
        hold on
        
        icol = 0;
        for iatt = 1:2
            for idrug = 1:2
                
                icol = icol+1;
                h(icol) = boundedline(times, ...
                    squeeze(mean(mean(samples(idx_unit,idrug,idx_attention==iatt,:),3),1)), ...
                    squeeze(std(mean(samples(idx_unit,idrug,idx_attention==iatt,:),3),[],1))/sqrt(length(find(idx_unit))), ...
                    'cmap', colors(iatt,idrug,:), 'alpha');
                set(h(icol), 'LineWidth', fSet.LineWidth)
                
                text_legend{icol} = sprintf('Attend %s, Drug %s', label_attention{iatt}, label_drug_onoff{idrug});
                
            end
        end
        xlim(xlim2use{ievent})
        ylim(ylim2use)
        
        if ievent==1
            ylabel('Normalized firing rate (%)', 'FontSize', fSet.Fontsize_text)
            
            text(200, 30, sprintf('n=%d', length(find(idx_unit))))
        elseif ievent==3
            if 0%irow==1
                hleg = legend(h, legend_text, 'FontSize', fSet.Fontsize_text, 'Box', 'Off');
                hleg.Position = hleg.Position + [0.05 0.15 0 0];
            end
        end
        
        if irow==nrow
            xlabel(xlabel_event{ievent}, 'FontSize', fSet.Fontsize_text)
        end
    end
    
    
end
hleg = legend(h, text_legend, 'FontSize', fSet.Fontsize_text, 'Box', 'Off', 'Location', 'NorthWest');
hleg.Position = hleg.Position + [-0.04 0.53 0 0];


% plot example cell, plot raster for example cell, only drug/attention
% ------------------------------------------------------------------------
% plot param
iplot = iplot+1;
axH = subtightplot(nrow, ncol, idx_subplot{iplot}, fSet.subplotGap, fSet.subplotMargin, fSet.subplotMargin);
hold on
plotj_initAx(fSet, 'axlabel', 2, 'axlabelDisplacement', [0.07, 0.00]);
yoffset = 100;


% load data
event2plot = 3;
% rec2plot = {'W','2015-10-16'};
rec2plot = {'W','2015-10-01'};
unit2plot = 1;

idx_unit = strcmpi(rec2plot(1), unitlist.Subject) ...
    & strcmpi(rec2plot(2), unitlist.Date) ...
    & (unitlist.unit==unit2plot);

loadfilename = fullfile(path_data, unitlist.Subject(idx_unit,:), unitlist.Date(idx_unit,:), 'unit.mat');
unit = load(loadfilename);
loadfilename = fullfile(path_data, unitlist.Subject(idx_unit,:), unitlist.Date(idx_unit,:), 'trialdata.mat');
load(loadfilename);
[trialdata, unit] = remove_excluded_trials(trialdata, unit, unit2plot);

% select only trial window for which this unit has spikes
sps = unit.Dim1Align;

% raster for each condition
idx_attend = idx_attention([trialdata.cond_num])';
idx_drug = [trialdata.drug]' + 1;
ymax = 0;
for iatt = 1:length(unique(idx_attention))
    for idrug = 1:2
        
        % get trial indices
        trial_index = (...
            idx_drug==idrug) ...
            & idx_attend==iatt;
        
        [xspikes, yspikes] = spike_raster( sps(:,trial_index), xlim2use{event2plot}, unit2plot);
        
        x = rate_PSTH.time;
        y = squeeze(mean(rate_PSTH.samples(idx_unit,idrug,idx_attention==iatt,:),3));
        h(icol) = plot(x, y, 'Color', colors(iatt,idrug,:), 'linew', fSet.LineWidth);
        
        plot(xspikes, yspikes+yoffset+ymax, 'Color', colors(iatt,idrug,:), 'linew', 1.5)
        
        ymax = ymax + max(yspikes(:));
    end
end
plot([-500 0], [10 10], 'Color', [0.5 0.5 0.5], 'linew', 5)
set(gca,'YTick',[0 yoffset],'YTicklabel',[0 round(rate_PSTH.maxhist(idx_unit)*1000)]);
xlabel('Time from first-dimming (ms)', 'FontSize', fSet.Fontsize_text)
ylabel('Firing rate (Hz)', 'FontSize', fSet.Fontsize_text)

xlim([-800 200])
ylim2use = yoffset+ymax+80;
ylim([0 ylim2use])

p = [rate_ANOVA.selectivity.att(idx_unit) rate_ANOVA.selectivity.dru(idx_unit) rate_ANOVA.selectivity.("att*dru")(idx_unit)];
[pstring,starstring] = get_significance_strings(p, 'rounded', 0, 'factorstring', {'attention', 'drug', 'interaction'});
text(axH, -750, ylim2use*0.95, pstring, 'FontSize', fSet.Fontsize_text)

% save
savefigname = fullfile(path_population, sprintf('Population_hist'));
plotj_saveFig(savefigname, {'png', 'svg', 'pdf'})



%% plot histogram cell-widths

binsize = 25;
binedges = 50:binsize:600;

assert(any(ismember(binedges,waveform_cutoff)), 'threshold needs to be included in binedges')

[fH, fSet] = plotj_initFig('width', 20, 'height', 12, 'Journal',plot_conventions);
fSet.subplotGap = fSet.subplotGap .* [1 0.8];

ncol=4;
nrow=2;

subtightplot(nrow, ncol, 1, fSet.subplotGap, fSet.subplotMargin, fSet.subplotMargin)
plotj_initAx(fSet, 'axlabel', 1, 'axlabelDisplacement', [0.07, 0.01]);
hold on

x = waveform.time-(7*6*5.4);
plot(x, waveform.waveform(unit_class==1, :)*-1, 'Color', [fSet.colors(1,:) 0.5])
plot(x, waveform.waveform(unit_class==2, :)*-1, 'Color', [fSet.colors(2,:) 0.5])

xlabel(['Time from peak of spike (' plotj_symbol('mu') 's)'], 'FontSize', fSet.Fontsize_text)
ylabel(['Normalized voltage'], 'FontSize', fSet.Fontsize_text)

xlim(minmax(waveform.time) + [0 -200] - 7*6*5.4 )
axis square

subtightplot(nrow, ncol, 2, fSet.subplotGap, fSet.subplotMargin, fSet.subplotMargin)
plotj_initAx(fSet, 'axlabel', 2, 'axlabelDisplacement', [0.07, 0.01]);
hold on

h = plotj_hist({waveform.peak_to_trough_time(unit_class==1), waveform.peak_to_trough_time(unit_class==2)}, ...
    'bins', binedges, 'histstyle', 'stairs',...
    'plotmean', 0, 'LineWidth', fSet.LineWidth);

xlabel(['Peak to through time (' char(181) 's)'], 'FontSize', fSet.Fontsize_text)
ylabel('Number of cells', 'FontSize', fSet.Fontsize_text)
axis square

hleg = legend(h, {'Narrow','Broad'}, 'FontSize', fSet.Fontsize_text, 'Box', 'Off', 'AutoUpdate', 'Off', 'Location', 'East');
hleg.Position = hleg.Position + [0.08 0 0 0];

% do hartigan's dip test to see whether there is a significant dip in this
% distribution
% [dip, p_value, xlow,xup]=HartigansDipSignifTest(allSpikeWaveform.width(idx_unit),5000);
[dip, p_value, xlow, xup, boot_dip]=CalibratedHartigansDipSignifTest(waveform.peak_to_trough_time, 10000);


p_string = get_significance_strings(p_value, 'round', 0);
h_text(1) = text(400, 13, {'Calibrated Hartigan''s dip test:', p_string}, 'FontSize', fSet.Fontsize_text);

savefigname = fullfile(path_population, sprintf('Population_broadNarrow'));
plotj_saveFig(savefigname, {'png', 'svg'})

%% mean rate, Fano factor, drugs

% selectivity_criterium = 'none';
selectivity_criterium = 'dru';
% selectivity_criterium = 'att&dru';

ncol = 2 * (length(label_drug));
nrow = 2; % rate, FF


idx_subplot = [1 2 ncol+1 ncol+2 ; 3 4 ncol+3 ncol+4] ;
idx_axlabel = [1 NaN 2 NaN ; 3 NaN 4 NaN]+2;

[fH, fSet] = plotj_initFig('figNum', 2, 'width', 20, 'height', 12, 'Journal',plot_conventions);
fSet.subplotGap = fSet.subplotGap .* [1 0.8];

[P, h_text] = deal( NaN(length(label_drug), 2, 2, 2) ); % drugtype, actvity-type, attention, drug-offon,
for idrug = 1:length(label_drug)
    
    % unit selection
    idx_unit = get_unit_selectivity(unitlist, rate_ANOVA.selectivity, selectivity_criterium, {'drug',label_drug(idrug)});
    
    iplot=0;
    for irow = 1:nrow
        
        if irow == 1
            % mean rate
            datatype = 'rate';
            data2plot = rate_summary.rate;
            ylim2use = [0.7 100];
            
            scale ='Log';
            x_text = 0.005;
            y_text = [0.9 0.4];
        elseif irow == 2
            % Fano factor
            datatype = 'FF';
            data2plot = rate_summary.FF;
            ylim2use = [0 10];
            
            scale ='Linear';
            x_text = 0.1;
            y_text = [0.95 0.1];
        end
        
        for iatt = 1:2
            iplot = iplot+1;
            
            subtightplot(nrow, ncol, idx_subplot(idrug, iplot), fSet.subplotGap, fSet.subplotMargin, fSet.subplotMargin)
            plotj_initAx(fSet, 'axlabel', idx_axlabel(idrug, iplot), 'axlabelDisplacement', [0.05, 0.02]);
            hold on
            
            title(sprintf('Attend %s', label_attention{iatt}), 'FontSize', fSet.Fontsize_title)
            
            for unittype = 1:2
                unit2plot = idx_unit & unit_class==unittype;
                
                % extract data
                tmp_data = squeeze(mean(data2plot(unit2plot,:,idx_attention==iatt),3));
                
                plotj_scatter(tmp_data, ...
                    'markerStyle', {'o'}, 'MarkerSize', markersize, ...
                    'markerFaceColor', fSet.colors(unittype,:), 'markerFaceAlpha', 0.5, ...
                    'markerEdgeColor', fSet.colors(unittype,:), 'markerEdgeAlpha', 0.5, ...
                    'axislimit', ylim2use);
                
                % stats
                P(idrug,irow,iatt,unittype) = compare_means(tmp_data(:,1), tmp_data(:,2), 1, 'rank');
                
                % effect size
                d = computeCohen_d(tmp_data(:,2), tmp_data(:,1), 'paired');
                
                % mean difference
                delta_data = tmp_data(:,2)-tmp_data(:,1);
                
                axis square
                p_string = get_significance_strings(P(idrug,irow,iatt,unittype), 'rounded', 0);
                
                % print result
                fprintf('%s: %s, %s, attend %s: delta-%s %1.2f +- %1.2f, %s, Cohens''s d = %1.2f\n', ...
                    datatype, label_drug{idrug}, label_celltype{unittype}, label_attention{iatt}, datatype, mean(delta_data), std(delta_data)/sqrt(length(find(unit2plot))), p_string, d)
                
                % plot text
                tmp_x = get(gca,'xlim');
                tmp_y = get(gca,'ylim');
                
                x_pos = get_value_range(tmp_x, x_text);
                y_pos = get_value_range(tmp_y, y_text(1)-y_text(2)*(unittype-1));
                
                h_text(idrug,irow,iatt,unittype) = text(x_pos, y_pos, ...
                    sprintf('%s', p_string), ...
                    'Color', fSet.colors(unittype,:));
                
                if iatt==2
                    
                    x_pos = get_value_range(tmp_x, 0.6);
                    y_pos = get_value_range(tmp_y, 0.1+0.1*(unittype-1));
                    
                    text(x_pos, y_pos, sprintf('%s (n=%d)', label_celltype{unittype}, length(find(unit2plot))), 'Color', fSet.colors(unittype,:));
                end
            end
            
            
            xlim(ylim2use)
            ylim(ylim2use)
            set(gca,'YScale', scale, 'XScale', scale)
            % %     set(gca,'', [1 10 100])
            
            if irow == 1
                set(gca,'Xtick', [1 10 100], 'XTickLabel', [1 10 100])
                set(gca,'Ytick', [1 10 100], 'YTickLabel', [1 10 100])
                
                xlabel('Firing rate no drug (Hz)','FontSize', fSet.Fontsize_text)
                
                if idrug==1 && iatt==1
                    ylabel('Firing rate drug (Hz)','FontSize', fSet.Fontsize_text)
                end
            elseif irow == 2
                xlabel('Fano Factor no drug','FontSize', fSet.Fontsize_text)
                if idrug==1 && iatt==1
                    ylabel('Fano Factor drug','FontSize', fSet.Fontsize_text)
                end
            end
        end
    end
end

[p_fdr, p_masked] = FDR(P, 0.05);
plotj_text_emphasise(h_text, p_masked, 'italic');
plotj_text_emphasise(h_text, p_masked, 'bold');


savefigname = fullfile(path_population, sprintf('mean_rate_FF_drug'));
plotj_saveFig(savefigname, {'png', 'svg'})


%% Attend ROC, drug/no drug, narrow/broad. DrugMI-ejecCurrent
datatype = 'MI';


selectivity_criterium = 'none';
selectivity_criterium = 'att&dru';
% selectivity_criterium = 'dru';
% selectivity_criterium = 'dru';
% selectivity_criterium = 'att';

ncol = length(label_drug);
nrow = 2;

[fH, fSet] = plotj_initFig('width', 20, 'height', 15, 'Journal',plot_conventions);
fSet.subplotGap = fSet.subplotGap.*[1.2 .8];

axislim = [0.25 0.9];
[P, h_text] = deal( NaN(length(label_drug), 2) ); % drugtype, unit-type,
[P_roc, h_text_roc] = deal( NaN(length(label_drug), 2) ); % drugtype, unit-type,
iplot=0;
for idrug = 1:length(label_drug)
    iplot = iplot+1;
    
    % unit selection
    idx_unit = get_unit_selectivity(unitlist, rate_ANOVA.selectivity, selectivity_criterium, {'drug',label_drug(idrug)});
    idx_unit_interaction = get_unit_selectivity(unitlist, rate_ANOVA.selectivity, 'att*dru', {'drug',label_drug(idrug)});
    
    % subplot axis
    h_ax_subplot = subtightplot(nrow, ncol, iplot, fSet.subplotGap, fSet.subplotMargin, fSet.subplotMargin);
    plotj_initAx(fSet, 'axlabel', idrug, 'axlabelDisplacement', [0.01, 0.02]);
    hold on
    
    title([label_drug{idrug} label_drug_ext{idrug}], 'FontSize', fSet.Fontsize_title)
    
    % inset axis
    pos = get(h_ax_subplot, 'Position');
    h_ax_inset = axes('Position', [pos(1)+0.28 pos(2)+0.035 0.1 0.11]) ; % inset
    plotj_initAx(fSet);
    hold on

    
    for unittype = 1:2
        
        set(fH, 'currentaxes', h_ax_subplot);
        
        unit2plot = idx_unit & unit_class==unittype;
        unit2plot_interaction = idx_unit_interaction & unit_class==unittype;
        
        tmp_data = squeeze(mean(rate_ROC.roc_attend(unit2plot,:,idx_attention_roc==1),3));
        idx = tmp_data(:,1)<0.5;
        tmp_data(idx,:) = 1-tmp_data(idx,:);

        
        idx_interaction = unit2plot_interaction(unit2plot)+1;
        
        plotj_scatter(tmp_data, ...
            'markerStyle', {'o','v'}, 'MarkerSize', markersize, ...
            'markerFaceColor', fSet.colors(unittype,:), 'markerFaceAlpha', 0.5, ...
            'markerEdgeColor', fSet.colors(unittype,:), 'markerEdgeAlpha', 0.5, ...
            'axislimit', axislim);
        
        % stats
        P(idrug,unittype) = compare_means(tmp_data(:,1),tmp_data(:,2), 1, 'rank');
        
        axis square
        p_string = get_significance_strings(P(idrug,unittype), 'rounded', 0);
        
        % set text
        tmp_x = get(gca,'xlim');
        tmp_y = get(gca,'ylim');
        
        x_pos = get_value_range(tmp_x, 0.05);
        y_pos = get_value_range(tmp_y, 0.95 - 0.07 * (unittype-1));
        
        h_text(idrug,unittype) = text(x_pos, y_pos, ...
            sprintf('%s, %s (n=%d)', p_string, label_celltype{unittype}, length(find(unit2plot))), ...
            'Color', fSet.colors(unittype,:), 'FontSize', fSet.Fontsize_ax);
        
        %%% inset
        % att AUROC relative to 0.5
%         tmp_data = abs(tmp_data-0.5);
%         P_roc(idrug,unittype) = compare_means(tmp_data(:,1),tmp_data(:,2), 1, 'rank');

        P_roc(idrug,unittype) = compare_means(tmp_data(:,1),tmp_data(:,2), 1, 'rank');
% 
        bar_data = diff(tmp_data,1,2);
        P_roc(idrug,unittype) = compare_means(bar_data,0, 1, 'rank');
% 
        
        bar_y = mean(bar_data);
        bar_y_err = std(bar_data)/sqrt(length(bar_data));
        
        % plot
        set(fH, 'currentaxes', h_ax_inset);
        h_scatter = scatter(unittype + (rand(length(bar_data),1)-0.5)*0.1, bar_data, ...
            '.', 'MarkerFaceColor', fSet.colors(unittype,:),'MarkerEdgeColor', fSet.colors(unittype,:));
        h_scatter.MarkerFaceAlpha = 0.5;
        h_scatter.MarkerEdgeAlpha = 0.5;

        plot(unittype + [-0.25 0.25], [bar_y bar_y], 'Color', fSet.colors(unittype,:), 'linew', fSet.LineWidth+1)
        plot(unittype + [0 0], bar_y + [-bar_y_err bar_y_err], 'Color', fSet.colors(unittype,:), 'linew', fSet.LineWidth)
        
        p_string = get_significance_strings(P_roc(idrug,unittype), 'rounded', 0);

        % set text
        % tmp_y = get(gca,'ylim');
        %y_pos = get_value_range(tmp_y, 0.95);
        y_pos = 0.1*idrug;
        
        h_text_roc(idrug,unittype) = text(unittype+0.1, y_pos, ...
            p_string, ...
            'Color', fSet.colors(unittype,:), ...
            'FontSize', fSet.Fontsize_ax);
    end
    
    set(h_ax_inset, 'XTick', [1 2], 'XTickLabel', label_celltype)
    ylabel(h_ax_inset, [plotj_symbol('Delta') ' | AUROC-0.5 |'],'FontSize', fSet.Fontsize_text)

    xlabel(h_ax_subplot, 'Attention AUROC no drug','FontSize', fSet.Fontsize_text)
    if idrug==1
        ylabel(h_ax_subplot, 'Attention AUROC drug','FontSize', fSet.Fontsize_text)
    end
end
[p_fdr, p_masked] = FDR(P, 0.05);
[p_fdr, p_masked_corr] = FDR(P_roc, 0.05);

plotj_text_emphasise(h_text, p_masked, 'italic');
plotj_text_emphasise(h_text, p_masked, 'bold');







%%% drugMI/attAUROC against ejection current

switch datatype
    case 'MI'
        %         data2plot = rate_ROC.mi_drug(:,1);
        data2plot = mean(rate_ROC.mi_drug,2);
        
        unittype=[];
        
        ylabel2use = 'Drug modulation index';
        ylim2use = [];
    case 'AUC'
        data2plot = squeeze(diff(rate_ROC.roc_attend,1,2));
        data2plot = squeeze(mean(data2plot(:,idx_attention_roc==1),2)); % average over relevant auroc conditions
        
        unittype = [];
        
        ylabel2use = [plotj_symbol('Delta') ' attention AUROC'];
        ylim2use = [-0.25 0.3];
end

clear ax h
for idrug = 1:ncol
    clear text_legend
    iplot=iplot+1;
    
    idx_unit = get_unit_selectivity(unitlist, rate_ANOVA.selectivity, selectivity_criterium, {'drug',label_drug(idrug)});
    
    ax(idrug) = subtightplot(nrow, ncol, iplot, fSet.subplotGap, fSet.subplotMargin, fSet.subplotMargin);
    plotj_initAx(fSet, 'axlabel', iplot, 'axlabelDisplacement', [0.01, 0.02]);
    hold on
    
    title([label_drug{idrug} label_drug_ext{idrug}], 'FontSize', fSet.Fontsize_title)
    
    if ~isempty(unittype)
        unit2plot = idx_unit & unit_class==unittype;
    else
        unit2plot = idx_unit;
    end
    length(find(unit2plot))
    
    % ejectionCurrent
    ejection_current = nanmean(unitlist.EjectCurrent(unit2plot),2);
    xlabel_event = 'Ejection Current';
    xlim2use = [17 90];
    
    if length(unique(unit_class(unit2plot)))==1
        tmp_data = [ejection_current  data2plot(unit2plot)];
        hscat = plotj_scatter(tmp_data, ...
            'markerStyle', {'.'}, 'markerSize', markersize, ...
            'markerFaceColor', fSet.colors(unittype,:), 'markerFaceAlpha', 0.5, ...
            'axislimit', [], ...
            'unityLine', 0);
        text_legend{1} = [label_celltype{unittype} ' (n=' num2str(length(find(unit_class(idx_unit)==unittype))) ')'];
    else
        tmp_data = [ejection_current  data2plot(unit2plot)];
        hscat = plotj_scatter(tmp_data, ...
            'dataIndex', unit_class(unit2plot), ...
            'markerStyle', {'o','o'}, 'MarkerSize', markersize, ...
            'markerFaceColor', fSet.colors(1:2,:), 'markerFaceAlpha', [0.5 0.5], ...
            'markerEdgeColor', fSet.colors(1:2,:), 'markerEdgeAlpha', [0.5 0.5], ...
            'axislimit', [], ...
            'unityLine', 0);
        text_legend{1} = [label_celltype{1} ' (n=' num2str(length(find(unit_class(idx_unit)==1))) ')'];
        text_legend{2} = [label_celltype{2} ' (n=' num2str(length(find(unit_class(idx_unit)==2))) ')'];
    end
    
    data_table = table(ejection_current, data2plot(unit2plot), 'VariableNames', {'x', 'y'});
    [stats2report, predict_data, stats, model] = fitlme_singleVar_sequential(data_table, 0);
    
    switch datatype
        case 'MI'
            filename = fullfile(path_population, sprintf('doseResponse_drugMI_%s_%s.csv', label_drug{idrug}, selectivity_criterium));
        case 'AUC'
            filename = fullfile(path_population, sprintf('doseResponse_attAUROC_%s_%s.csv', label_drug{idrug}, selectivity_criterium));
    end
    writetable(data_table, filename)
    
    fprintf('Chi(%d) %1.2f, p = %1.3f\n', stats2report.deltaDF, stats2report.LRStat, stats2report.pValue)
    
    if any(strcmp('fit_mean', predict_data.Properties.VariableNames))
        h = plot(predict_data.x, predict_data.fit_mean, 'linewidth', fSet.LineWidth, 'color', [0.3 0.3 0.3]);
        h = plot(predict_data.x, predict_data.fit_CI_low, '--', 'linewidth', fSet.LineWidth_in, 'Color', [0.3 0.3 0.3]);
        h = plot(predict_data.x, predict_data.fit_CI_high, '--', 'linewidth', fSet.LineWidth_in, 'Color', [0.3 0.3 0.3]);
        [pstring_chi,starstring] = get_significance_strings(stats2report.pValue, 'rounded', 0);
        
        cond1 = exist([filename(1:end-4) '_R.csv'], 'file');
        cond2 = strcmpi(model,'quadratic');
        if cond1 && cond2
            utable = readtable([filename(1:end-4) '_R.csv']);
            
            if utable.usig==1
                pstring_chi = [pstring_chi ' (U^{+})'];
            else
                pstring_chi = [pstring_chi ' (U^{-})'];
            end
        end
        
    else
        pstring_chi = '';
    end
    
    
    switch model
        case 'linear'
            betaString = ['\beta' '_{1} = ' num2str(stats.Estimate, '%1.2e')];
            %         model_color = fSet.colors(2,:);
        case 'quadratic'
            betaString = ['\beta' '_{2} = ' num2str(stats.Estimate, '%1.2e')];
        otherwise
            betaString = [];
    end
    
    axis square
    xlabel(xlabel_event, 'FontSize', fSet.Fontsize_text)
    if idrug==1
        ylabel(ylabel2use, 'FontSize', fSet.Fontsize_text)
    end
    
    if ~isempty(ylim2use)
        ylim(ylim2use)
    end
    if ~isempty(xlim2use)
        xlim(xlim2use)
    end
    
    % plot text indicating number of units
    set(gca, 'Units', 'normalized');
    
    tmp_x = get(gca,'xlim');
    tmp_y = get(gca,'ylim');
    
    x_pos = get_value_range(tmp_x, 0.6);
    y_pos(1) = get_value_range(tmp_y, 0.8);
    y_pos(2) = get_value_range(tmp_y, 0.9);
    
    text(x_pos, y_pos(1), text_legend{1}, 'FontSize', fSet.Fontsize_text, 'Color', fSet.colors(1,:))
    text(x_pos, y_pos(2), text_legend{2}, 'FontSize', fSet.Fontsize_text, 'Color', fSet.colors(2,:))

    % set ax param
    set(gca,'YDir','reverse')
    
    tmp_x = get(gca,'xlim');
    tmp_y = get(gca,'ylim');
    
    x_pos = get_value_range(tmp_x, 0.05);
    y_pos = get_value_range(tmp_y, 0.1);
    
    text(x_pos, y_pos, {betaString, pstring_chi }, 'fontsize', fSet.Fontsize_text, 'color', [0.3 0.3 0.3])
    
end

savefigname = fullfile(path_population, sprintf('attAuroc_DrugResponseCurve'));
plotj_saveFig(savefigname, {'png', 'svg'})


%% DrugMI-ejecCurrent across selectivity indices

datatype = 'MI';
% datatype = 'AUC';


selectivity_criteria = {'none', 'att', 'dru', 'att&dru'};
label_criteria = {'No subselection', 'Attention-selective', 'Drug-selective', {'Attention & Drug', 'selective'}};

ncol = length(selectivity_criteria);
nrow = length(label_drug)+1;
iplot = 0;

[fH, fSet] = plotj_initFig('width', 20, 'height', 16, 'Journal',plot_conventions);
fSet.subplotGap = fSet.subplotGap.*[1 .75];

idx_subplot = [1:ncol ; (ncol+1):ncol*2; (ncol*2+1):ncol*3];
%%% drugMI/attAUROC against ejection current
[P,h_text] = deal(zeros(length(selectivity_criteria), length(label_drug)));
[P_population,h_text_population] = deal(zeros(length(selectivity_criteria),1));
for icrit = 1:length(selectivity_criteria)
    
    selectivity_criterium = selectivity_criteria{icrit};
    
    switch datatype
        case 'MI'
            %         data2plot = rate_ROC.mi_drug(:,1);
            data2plot = mean(rate_ROC.mi_drug,2);
            
            unittype=[];
            
            ylabel2use = 'Drug modulation index';
            ylim2use = [];
        case 'AUC'
            data2plot = squeeze(diff(rate_ROC.roc_attend,1,2));
            data2plot = squeeze(mean(data2plot(:,idx_attention_roc==1),2)); % average over relevant auroc conditions
            
            unittype = [];
            
            ylabel2use = [plotj_symbol('Delta') ' attention AUROC'];
            ylim2use = [-0.25 0.3];
    end
    
    clear ax h
    for idrug = 1:length(label_drug)
        clear text_legend
        iplot=iplot+1;
        
        idx_unit = get_unit_selectivity(unitlist, rate_ANOVA.selectivity, selectivity_criterium, {'drug',label_drug(idrug)});
        
        ax(idrug) = subtightplot(nrow, ncol, idx_subplot(iplot), fSet.subplotGap, fSet.subplotMargin, fSet.subplotMargin);
        if idrug==1
            plotj_initAx(fSet, 'axlabel', icrit);
        else
            plotj_initAx(fSet);
        end
        hold on
        
        if idrug==1
            % title([label_drug{idrug} label_drug_ext{idrug}], 'FontSize', fSet.Fontsize_title)
            h_title = title(label_criteria{icrit}, 'FontSize', fSet.Fontsize_title);
            set(h_title, 'Units', 'normalized');
            h_title.Position = h_title.Position+[0 0.2 0];
        end
        if ~isempty(unittype)
            unit2plot = idx_unit & unit_class==unittype;
        else
            unit2plot = idx_unit;
        end
        
        % ejectionCurrent
        ejection_current = nanmean(unitlist.EjectCurrent(unit2plot),2);
        xlabel_event = 'Ejection Current';
        xlim2use = [17 95];
        
        if length(unique(unit_class(unit2plot)))==1
            tmp_data = [ejection_current  data2plot(unit2plot)];
            hscat = plotj_scatter(tmp_data, ...
                'markerStyle', {'.'}, 'markerSize', markersize, ...
                'markerFaceColor', fSet.colors(unittype,:), 'markerFaceAlpha', 0.5, ...
                'axislimit', [], ...
                'unityLine', 0);
            text_legend{1} = [label_celltype{unittype} ' (n=' num2str(length(find(unit_class(idx_unit)==unittype))) ')'];
        else
            tmp_data = [ejection_current  data2plot(unit2plot)];
            hscat = plotj_scatter(tmp_data, ...
                'dataIndex', unit_class(unit2plot), ...
                'markerStyle', {'o','o'}, 'MarkerSize', markersize, ...
                'markerFaceColor', fSet.colors(1:2,:), 'markerFaceAlpha', [0.5 0.5], ...
                'markerEdgeColor', fSet.colors(1:2,:), 'markerEdgeAlpha', [0.5 0.5], ...
                'axislimit', [], ...
                'unityLine', 0);
            text_legend{1} = [label_celltype{1} ' (n=' num2str(length(find(unit_class(idx_unit)==1))) ')'];
            text_legend{2} = [label_celltype{2} ' (n=' num2str(length(find(unit_class(idx_unit)==2))) ')'];
        end
        
        data_table = table(ejection_current, data2plot(unit2plot), 'VariableNames', {'x', 'y'});
        [stats2report, predict_data, stats, model] = fitlme_singleVar_sequential(data_table, 0);
        P(icrit,idrug) = stats2report.pValue;
                
        switch datatype
            case 'MI'
                filename = fullfile(path_population, sprintf('doseResponse_drugMI_%s_%s.csv', label_drug{idrug}, selectivity_criteria{icrit}));
            case 'AUC'
                filename = fullfile(path_population, sprintf('doseResponse_attAUROC_%s_%s.csv', label_drug{idrug}, selectivity_criteria{icrit}));
        end
        writetable(data_table, filename);
        
        fprintf('Chi(%d) %1.2f, p = %1.3f\n', stats2report.deltaDF, stats2report.LRStat, stats2report.pValue)
        
        if any(strcmp('fit_mean', predict_data.Properties.VariableNames))
            h = plot(predict_data.x, predict_data.fit_mean, 'linewidth', fSet.LineWidth, 'color', [0.3 0.3 0.3]);
            h = plot(predict_data.x, predict_data.fit_CI_low, '--', 'linewidth', fSet.LineWidth_in, 'Color', [0.3 0.3 0.3]);
            h = plot(predict_data.x, predict_data.fit_CI_high, '--', 'linewidth', fSet.LineWidth_in, 'Color', [0.3 0.3 0.3]);
            [pstring_chi,starstring] = get_significance_strings(stats2report.pValue, 'rounded', 0);
            
            cond1 = exist([filename(1:end-4) '_R.csv'], 'file');
            cond2 = strcmpi(model,'quadratic');
            if cond1 && cond2
                utable = readtable([filename(1:end-4) '_R.csv']);
                
                if utable.usig==1
                    pstring_chi = [pstring_chi ' (U^{+})'];
                else
                    pstring_chi = [pstring_chi ' (U^{-})'];
                end
            end
            
        else
            pstring_chi = '';
        end
        
        switch model
            case 'linear'
                betaString = ['\beta' '_{1} = ' num2str(stats.Estimate, '%1.2e')];
                %         model_color = fSet.colors(2,:);
            case 'quadratic'
                betaString = ['\beta' '_{2} = ' num2str(stats.Estimate, '%1.2e')];
            otherwise
                betaString = [];
        end
        
        axis square
        if idrug==length(label_drug)
            xlabel(xlabel_event, 'FontSize', fSet.Fontsize_text)
        end
        if icrit==1
            ylabel({label_drug{idrug},ylabel2use}, 'FontSize', fSet.Fontsize_text)
        end
        
        if ~isempty(ylim2use)
            ylim(ylim2use)
        end
        if ~isempty(xlim2use)
            xlim(xlim2use)
        end
        
        % plot text indicating number of units
        set(gca, 'Units', 'normalized');
        
        tmp_x = get(gca,'xlim');
        tmp_y = get(gca,'ylim');
        
        x_pos = get_value_range(tmp_x, 0.6);
        y_pos(1) = get_value_range(tmp_y, 0.8);
        y_pos(2) = get_value_range(tmp_y, 0.9);
        
        text(x_pos, y_pos(1), text_legend{1}, 'FontSize', fSet.Fontsize_text, 'Color', fSet.colors(1,:))
        text(x_pos, y_pos(2), text_legend{2}, 'FontSize', fSet.Fontsize_text, 'Color', fSet.colors(2,:))
        
        % set ax param
        set(gca,'YDir','reverse')
        
        % plot text with stats
        tmp_x = get(gca,'xlim');
        tmp_y = get(gca,'ylim');
        
        switch label_drug{idrug}
            case 'Dopamine'
                x_pos = get_value_range(tmp_x, 0.65);
                y_pos = get_value_range(tmp_y, 0.05);
            case 'SCH23390'
                x_pos = get_value_range(tmp_x, 0.05);
                y_pos = get_value_range(tmp_y, 0);
        end
        
        h_text(icrit,idrug) = text(x_pos, y_pos, {betaString, pstring_chi}, 'fontsize', fSet.Fontsize_text, 'color', [0.3 0.3 0.3]);
        
    end
    
    
    
    %%% double fit of lme
    ejection_current = nanmean(unitlist.EjectCurrent,2);
    idx_unit = get_unit_selectivity(unitlist, rate_ANOVA.selectivity, selectivity_criterium);

    tmp_unitlist = unitlist(idx_unit,:);
    idx_drug = strcmpi(tmp_unitlist.Drug, label_drug(2))+1;
    
    data_table = table(ejection_current(idx_unit), data2plot(idx_unit), idx_drug, 'VariableNames', {'x', 'y', 'idx'});

    switch datatype
        case 'MI'
            filename = fullfile(path_population, sprintf('doseResponse_drugMI_%s.csv', selectivity_criteria{icrit}));
        case 'AUC'
            filename = fullfile(path_population, sprintf('doseResponse_attAUROC_%s.csv', selectivity_criteria{icrit}));
    end
    writetable(data_table, filename);
    
    
    data_table.polyx1 = data_table.x;
    data_table.polyx2 = data_table.x.^2;
    
    data_table.mixed_poly = data_table.polyx1;
    data_table.mixed_poly(idx_drug==1) = data_table.polyx2(idx_drug==1);
    
    
    lme1 = fitlme(data_table,['y ~ 1']);% fit constant
    lme2 = fitlme(data_table,['y ~ polyx1']); % fit linear model
    lme3 = fitlme(data_table,['y ~ polyx1 + polyx2']); % fit quadratic model (this includes the linear model), note this doesn't have orthogonal contrasts!!    
    lme4 = fitlme(data_table,['y ~ polyx1 + polyx2 + mixed_poly']); % fit mixed linear/quadratic coefficients
    
    % compare the model fits
    comp1 = compare(lme1, lme2, 'CheckNesting',true); % compare linear to constant
    comp2 = compare(lme2, lme3, 'CheckNesting',true); % compare linear to constant
    comp3 = compare(lme3, lme4, 'CheckNesting',true); % compare linear to constant
    
    P_population(icrit) = comp3.pValue(2);
    
    % print results
    fprintf('\nPopulation: \n')
    fprintf('Linear: Chi(%d) %1.2f, p = %1.3f\n', 1, comp1.LRStat(2), comp1.pValue(2))
    fprintf('Quadratic: Chi(%d) %1.2f, p = %1.3f\n', 1, comp2.LRStat(2), comp2.pValue(2))
    fprintf('Mixed: Chi(%d) %1.2f, p = %1.3f\n', 1, comp3.LRStat(2), comp3.pValue(2))

    % plot text
    iplot=iplot+1;    
    subtightplot(nrow, ncol, idx_subplot(iplot), fSet.subplotGap, fSet.subplotMargin, fSet.subplotMargin);
    plotj_initAx(fSet);
    hold on
    axis off
    [pstring_chi] = get_significance_strings(comp3.pValue(2), 'rounded', 0);
    [pstring_coef1] = get_significance_strings(lme4.Coefficients.pValue(2), 'rounded', 0);
    [pstring_coef2] = get_significance_strings(lme4.Coefficients.pValue(3), 'rounded', 0);
    [pstring_coef3] = get_significance_strings(lme4.Coefficients.pValue(4), 'rounded', 0);

    
    
    text2plot = { ...
        sprintf('%s(1) = %1.2f, %s', '\chi', comp3.LRStat(2), pstring_chi), ... # final model fit
        '', ...
        sprintf('%s_{1} = %1.2e, %s', '\beta', lme4.Coefficients.Estimate(2), pstring_coef1), ... # first coefficient
        sprintf('%s_{2} = %1.2e, %s', '\beta', lme4.Coefficients.Estimate(3), pstring_coef2), ... # first coefficient
        sprintf('%s_{3} = %1.2e, %s', '\beta', lme4.Coefficients.Estimate(4), pstring_coef3), ... # first coefficient
        };
    %

    h_text_population(icrit,1) = text(0,0,text2plot,'Interpreter','tex');
    ylim([-0.5 0.5])
    
    if icrit==1
        text(-0.2, 0, ['"Population Model"'], 'Rotation', 90, 'HorizontalAlignment', 'center')
    end

end

% multiple comparison correction, individual models
[p_fdr, p_masked] = FDR(P, 0.05);
plotj_text_emphasise(h_text, p_masked, 'italic', 2);
plotj_text_emphasise(h_text, p_masked, 'bold', 2);

% multiple comparison correction, population models
[p_fdr, p_masked] = FDR(P_population, 0.05);
plotj_text_emphasise(h_text_population, p_masked, 'italic', 1);
plotj_text_emphasise(h_text_population, p_masked, 'bold', 1);


% save
savefigname = fullfile(path_population, sprintf('DrugResponseCurve_selection'));
plotj_saveFig(savefigname, {'png', 'svg'})

%% plot behaviour

mfactor = 1000; %multiplication_factor
colors = get_colors('att_drug');

ncol = length(label_drug);
nrow = 1;

ylim2plot = [250 550];

[fH, fSet] = plotj_initFig('width', 25, 'height', 10, 'Journal',plot_conventions);

clear legendText h
for idrug = 1:ncol
    
    idx_rec = strcmpi(recordinglist.Drug, label_drug(idrug));
    %
    subtightplot(nrow, ncol, idrug, fSet.subplotGap, fSet.subplotMargin, fSet.subplotMargin)
    plotj_initAx(fSet, 'axlabel', idrug, 'axlabelDisplacement', [0.07, 0.01]);
    hold on
    
    title([label_drug{idrug} label_drug_ext{idrug}], 'FontSize', fSet.Fontsize_title)
    
    %%% plot errorbar
    iattDrug = 0;
    RT_table = [];
    for iatt = 1:2
        for idrug_in = 1:2
            iattDrug = iattDrug+1;
            
            x = iattDrug;
            xjitter = (rand(length(find(idx_rec)),1)-0.5)/3;
            y = squeeze(mean(RT.RT(idx_rec,idrug_in,idx_attention==iatt),3)) * mfactor;
            
            plot(x + xjitter, y, '.', ...
                'Color', [0.2 0.2 0.2 0.5], ...
                'MarkerSize', 10)
            
            h(iattDrug) = plotj_errorBar(y, 'x2plot', iattDrug);
            h(iattDrug).FaceColor = colors(iatt,idrug_in,:);
            h(iattDrug).FaceAlpha = 0.7;
            
            text_legend{iattDrug} = sprintf('Attend %s, Drug %s', label_attention{iatt}, label_drug{idrug_in});
            
            %             plot(iattDrug + xjitter, squeeze(allRTc2plot(:,iatt,idrug_in)), '.', 'Color', attDrugColors(iatt,idrug_in,:))
            
            n = length(find(idx_rec));
            rec = (1:n)';
            cond = repmat(iatt,n,1);
            drug = repmat(idrug_in,n,1);
            RT_table = [RT_table ; table(rec,y/mfactor,cond,drug,'VariableNames',{'recording','RT','attention','drug'})];
        end
    end
    
    %%% do stats
    % check with linear mixed effect model
    lme = fitlme(RT_table,['RT ~ 1 + (1|recording)']);% fit constant
    lme2 = fitlme(RT_table,['RT ~ attention + (1|recording)']); % fit linear model with attention condition
    lme3 = fitlme(RT_table,['RT ~ attention + drug + (1|recording)']); %  fit linear model with states
    lme4 = fitlme(RT_table,['RT ~ attention + drug + attention:drug + (1|recording)']); % fit interaction
    
    [BETA,BETANAMES,STATS] = fixedEffects(lme4);
    
    % compare the model fits
    comp1 = compare(lme, lme2, 'CheckNesting',true); % compare linear to constant
    comp2 = compare(lme2, lme3, 'CheckNesting',true); %
    comp3 = compare(lme3, lme4, 'CheckNesting',true); %
    
    [p_fdr, p_masked] = FDR( [comp1.pValue(2) comp2.pValue(2) comp3.pValue(2)], 0.05);
    
    pstring = get_significance_strings([STATS.pValue(2:end)], 'rounded', 0, 'factorstring', {'attention', 'drug', 'interaction'});
    
    pstring{4} = ['n = ' num2str(length(find(idx_rec)))];
    
    text(0.5, 500, pstring, 'FontSize', fSet.Fontsize_text)
    
    if idrug==1
        ylabel('Reaction time (ms)', 'FontSize', fSet.Fontsize_text)
    end
    
    set(gca, 'xtick', [1:4], 'xticklabel', text_legend, 'xticklabelRotation', -25)
    
    ylim(ylim2plot)
end


savefigname = fullfile(path_population, sprintf('RT'));
plotj_saveFig(savefigname, {'png', 'svg'})


%% pupil timecourse

events = {'stim', 'cue', 'dim'};
xlabel_event = {'Time from stimulus onset (ms)', 'Time from cue onset (ms)', 'Time from first-dimming (ms)'};
xlim2use = {[-200 1000], [-200 1000], [-1000 200]};
ylim2use = [-2.7 0.5];

ncol = length(events); % events
nrow = 2;%length(label_drug)-1;

timewinsize = 200;
timestep = 10;

[fH, fSet] = plotj_initFig('width', 20, 'height', 10, 'Journal',plot_conventions);
fSet.subplotGap = fSet.subplotGap ./ [1 2];
clear ax h legendText
iplot = 0;
for idrug = 1:nrow
    
    for ievent = 1:length(events)
        iplot = iplot+1;
        
        
        [pupil_timeseries, reclist] = get_population_data(recordinglist, 'pupil_timeseries', path_data, events{ievent});
        
        if nrow==length(label_drug)
            idx_rec = strcmpi(reclist.Drug, label_drug(idrug));
        else
            idx_rec = true(height(reclist),1);
        end
        nrecs = length(find(idx_rec));
        
        pupil2plot = squeeze(nanmean(pupil_timeseries.pupil(idx_rec,:,:,:),3));
        time = pupil_timeseries.time*1000;
        
        timecenters = time(1):timestep:time(end);
        p_drug = NaN(1,length(timecenters));
        
        clear legendText
        
        h_ax = subtightplot(nrow, ncol, iplot, fSet.subplotGap, fSet.subplotMargin, fSet.subplotMargin);
        if ievent==1
            plotj_initAx(fSet, 'axlabel', idrug, 'axlabelDisplacement', [0.07, 0.01]);
        else
            plotj_initAx(fSet);
        end
        hold on
        
        if (ievent~=1)
            axh = get(gca);
            axh.YTick = [];
            axh.YAxis.Visible = 'off'; % remove y-axis
        end
        
        ilegend = 0;
        
        for idrug_in = 1:2
            ilegend = ilegend+1;
            
            y = squeeze(nanmean(pupil2plot(:,idrug_in,:),1));
            y_err = squeeze(nanstd(pupil2plot(:,idrug_in,:),[],1))/sqrt(size(pupil2plot,1));
            
            h(ilegend) = boundedline(...
                time, ...
                y, ...
                y_err, ...
                'alpha', 'cmap', fSet.colors(idrug_in,:));
            set(h(ilegend), 'LineWidth', fSet.LineWidth)
            
            text_legend{ilegend} = sprintf('Drug %s', label_drug_onoff{idrug_in});
            
        end
        
        for itc = 1:length(timecenters)
            
            timewin = timecenters(itc) + [-timewinsize/2 timewinsize/2];
            idx_time = dsearchn(time(:), timewin(:));
            
            tmp_data = squeeze(nanmean(pupil2plot(:,:,idx_time(1):idx_time(2)),3));
            
            p_drug(itc) = compare_means(tmp_data(:,1), tmp_data(:,2), 1, 'rank');
        end
        
        [p_fdr, p_masked] = FDR( p_drug, 0.05);
        plotj_brokenVector(timecenters, ylim2use(1)*0.95, p_masked, 'linewidth', 5, 'color', [0.5 0.5 0.5]);
        
        if iplot==1
            hleg = legend(h, text_legend, 'autoupdate', 'off', 'FontSize', fSet.Fontsize_text, 'Box', 'Off');
        end
        xlim(xlim2use{ievent})
        ylim(ylim2use)
        plot([0 0], ylim, 'k','linew',1)
        
        if ievent==1
            ylabel('Pupil diameter (z-score)', 'FontSize', fSet.Fontsize_text)
        end
        if idrug==nrow
            xlabel(xlabel_event{ievent}, 'FontSize', fSet.Fontsize_text)
        end
        
        if ievent~=1
            set(gca,'yticklabel',[])
        end
        if nrow>1
            if ievent==2
                title([label_drug{idrug} label_drug_ext{idrug}], 'FontSize', fSet.Fontsize_title)
            end
        end
        
        if ievent==1
            text(h_ax, 500, -0.5, sprintf('n=%d',nrecs), 'FontSize', fSet.Fontsize_text)
        end
        
    end
end
savefigname = fullfile(path_population, sprintf('Pupil_timecourse'));
plotj_saveFig(savefigname, {'png', 'svg'})

%% Pupil drug/no drug, scatter
%

colors = get_colors('drug');

events = {'baseline', 'stim', 'cue', 'dim'};
event_title = {'Baseline','Stimulus evoked','Cue evoked','First-dimming'};

ncol = length(events); % events
nrow = 1;

xtext = [22 22];
ytext = [-0.45 -0.32];


[fH, fSet] = plotj_initFig('width', 20, 'height', 10, 'Journal',plot_conventions);
fSet.subplotGap = fSet.subplotGap/1.5;
clear ax h legendText

[P, h_text] = deal(NaN(length(events), length(label_drug)));
for ievent = 1:length(events)
    
    [pupil, reclist] = get_population_data(recordinglist, 'pupil_windows', path_data, events{ievent});
    
    axlim = [];
    
    subtightplot(nrow, ncol, ievent, fSet.subplotGap, fSet.subplotMargin, fSet.subplotMargin)
    plotj_initAx(fSet, 'axlabel', ievent+2, 'axlabelDisplacement', [0.05, 0.13]);
    hold on
    title(event_title{ievent}, 'FontSize', fSet.Fontsize_title)
    
    % significance testing
    text_p = cell(length(label_drug),1);
    for idrug = 1:length(label_drug)
        
        idx_rec = strcmpi(reclist.Drug, label_drug(idrug));
        
        pupil2plot = squeeze(nanmean(pupil.pupil(idx_rec,:,:),3));
        
        P(ievent,idrug) = compare_means(pupil2plot(:,1), pupil2plot(:,2), 1, 'rank');
        
        text_p{idrug} = get_significance_strings(P(ievent,idrug), 'rounded', 0);
        
        text_legend{idrug} = sprintf('%s (n=%d)', label_drug{idrug}, length(find(idx_rec)));
    end
    
    
    pupil2plot = squeeze(nanmean(pupil.pupil,3));
    h = plotj_scatter(pupil2plot, ...
        'dataIndex', idx_rec+1, ...
        'Markersize', markersize, 'markerStyle', {'o','s'}, ...
        'MarkerEdgeColor', colors, 'MarkerFaceColor', colors, ...
        'MarkerEdgeAlpha', [0.5 0.5] ,'MarkerFaceAlpha', [0.5 0.5]);
    
    %     set(gca,'ydir','reverse')
    axis square
    tmp_x = get(gca,'xlim');
    tmp_y = get(gca,'ylim');
    
    x_pos = get_value_range(tmp_x, 0.6);
    
    for idrug = 1:length(label_drug)
        
        y_pos = get_value_range(tmp_y, 0.2-0.1*(idrug-1));
        
        h_text(ievent,idrug) = text(x_pos, y_pos, text_p{idrug}, 'FontSize', fSet.Fontsize_text, 'color', colors(idrug,:));
    end
    
    if ievent==1
        ylabel({'Drug (z-score)'}, 'FontSize', fSet.Fontsize_text)
    end
    xlabel({'No drug (z-score)'}, 'FontSize', fSet.Fontsize_text)
    
    if ievent==4
        hleg = legend(h, text_legend, 'FontSize', fSet.Fontsize_text, 'Box', 'Off', 'Location', 'NorthEast');
        hleg.Position = hleg.Position + [0.05 0.3 0 0];
    end
    
end

[p_fdr, p_masked] = FDR(P, 0.05);
plotj_text_emphasise(h_text, p_masked, 'italic');
plotj_text_emphasise(h_text, p_masked, 'bold');

savefigname = fullfile(path_population, sprintf('Pupil_scatter'));
plotj_saveFig(savefigname, {'png', 'svg'})

