function generateReportLangloc_v2(obj, reportName)
    % Check if this is a langloc experiment, if not, return early
    langlocExperiments = {'LangLocVisual', 'LangLoc', 'MITLangloc', 'LangLocAudio', 'LangLocAudio-2','analysis'};
    
    if ~ismember(obj.experiment, langlocExperiments)
        fprintf('Skipping non-langloc experiment: %s\n', obj.experiment);
        return;
    end
    
    % Create the crunched folder if it doesn't exist
    crunchedFolder = fullfile(obj.crunched_file_path);
    if ~exist(crunchedFolder, 'dir')
        mkdir(crunchedFolder);
    end

    % Create a new PDF file in the crunched folder
    pdfFileName = fullfile(crunchedFolder, strcat(reportName, '.pdf'));
    if exist(pdfFileName, 'file')
        delete(pdfFileName);
    end

    % Start the PDF document with proper title page
    import mlreportgen.report.*
    import mlreportgen.dom.*

    tempFolder = tempdir;
    
    % Create report object with title page
    rpt = Report(pdfFileName, 'pdf');
    
    % Add title page
    tp = TitlePage();
    tp.Title = 'Language Localization Experiment Report';
    tp.Subtitle = 'Neural Data Analysis';
    tp.Author = 'Kumar Duraivel';
    tp.PubDate = datestr(now, 'dd-mmm-yyyy');
    add(rpt, tp);

    % General Experiment Information
    add(rpt, Heading1('General Experiment Information'));
    
    if(sum(contains(obj.elec_ch_type, 'seeg')))
        if(isfield(obj.stats.sig_hg_channel,'pChan_bip'))
            infoTable = {
                'Subject Name', char(obj.subject);
                'Experiment Name', char(obj.experiment);
                'Total Trials Completed', sprintf('%d', size(obj.events_table,1));
                'Total Electrodes Implanted', sprintf('%d', length(obj.elec_ch_label));
                'Total SEEG Electrodes', sprintf('%d', sum(contains(obj.elec_ch_type, 'seeg')));
                'Total ECoG Grid Electrodes', sprintf('%d', sum(contains(obj.elec_ch_type, 'ecog_grd')));
                'Total ECoG Strip Electrodes', sprintf('%d', sum(contains(obj.elec_ch_type, 'ecog_strip')));
                'Unipolar Contacts', sprintf('%d', length(obj.elec_ch_valid));
                'Bipolar Contacts', sprintf('%d', length(obj.bip_ch_label));
                'Significant Unipolar Electrodes (High Gamma)', sprintf('%d', sum(obj.stats.sig_hg_channel.h_fdr_05));
                'Significant Bipolar Electrodes (High Gamma)', sprintf('%d', sum(obj.stats.sig_hg_channel.h_bip_fdr_05));
            };
        else
            infoTable = {
                'Subject Name', obj.subject;
                'Experiment Name', obj.experiment;
                'Total Trials Completed', num2str(size(obj.events_table,1));
                'Total Electrodes Implanted', num2str(length(obj.elec_ch_label));
                'Total SEEG Electrodes', num2str(sum(contains(obj.elec_ch_type, 'seeg')));
                'Total ECoG Grid Electrodes', num2str(sum(contains(obj.elec_ch_type, 'ecog_grd')));
                'Total ECoG Strip Electrodes', num2str(sum(contains(obj.elec_ch_type, 'ecog_strip')));
                'Unipolar Contacts', num2str(length(obj.elec_ch_valid));
                'Bipolar Contacts', num2str(length(obj.bip_ch_label));
                'Significant Unipolar Electrodes (High Gamma)', num2str(sum( obj.stats.sig_hg_channel.h_fdr_05));
            };
        end
    else
        infoTable = {
            'Subject Name', obj.subject;
            'Experiment Name', obj.experiment;
            'Total Trials Completed', num2str(size(obj.events_table,1));
            'Total Electrodes Implanted', num2str(length(obj.elec_ch_label));
            'Total SEEG Electrodes', num2str(sum(contains(obj.elec_ch_type, 'seeg')));
            'Total ECoG Grid Electrodes', num2str(sum(contains(obj.elec_ch_type, 'ecog_grd')));
            'Total ECoG Strip Electrodes', num2str(sum(contains(obj.elec_ch_type, 'ecog_strip')));
            'Unipolar Contacts', num2str(length(obj.elec_ch_valid));
            'Bipolar Contacts', num2str(length(obj.bip_ch_label));
            'Significant Unipolar Electrodes (High Gamma)', num2str(sum(cellfun(@(x) any(x.h_sig_05), obj.stats.time_series.pSigChan)));
        };
    end
    % Ensure all entries are properly formatted strings
    for i = 1:size(infoTable, 1)
        for j = 1:size(infoTable, 2)
            if isnumeric(infoTable{i,j})
                infoTable{i,j} = sprintf('%g', infoTable{i,j});
            elseif ~ischar(infoTable{i,j}) && ~isstring(infoTable{i,j})
                infoTable{i,j} = char(string(infoTable{i,j}));
            end
        end
    end
    tbl = Table(infoTable);
    tbl.Style = {Border('solid'), ColSep('solid'), RowSep('solid')};
    tbl.TableEntriesStyle = {HAlign('left')};
    add(rpt, tbl);

    % Data Extraction
    acc = [obj.events_table.accuracy];
    rt = [obj.events_table.RT];
    cond = [obj.condition];

    % Performance Metrics
    add(rpt, Heading2('Performance Metrics'));

    % Overall Accuracy
    accuracy_percentage = mean(acc)*100;
    add(rpt, Paragraph(['Overall accuracy: ' num2str(accuracy_percentage, '%.1f') '%']));

    % Reaction Time Analysis
    rt_correct = rt(acc == 1);
    cond_correct = cond(acc == 1);

    % Process only langloc experiments
    switch obj.experiment
        case {'LangLocVisual','LangLoc','MITLangloc'}
            mean_rt_all = mean(rt_correct, 'omitnan');
            mean_rt_sentence = mean(rt_correct(strcmp(cond_correct, 'sentence')), 'omitnan');
            mean_rt_nonword = mean(rt_correct(strcmp(cond_correct, 'nonword')), 'omitnan');
        
            % Results Table
            tableContent = {
                'Condition', 'Mean RT (ms)';
                'Overall', sprintf('%6.1f ± %.1f', mean_rt_all*1000, std(rt_correct)*1000);
                'Sentence (S)', sprintf('%6.1f ± %.1f', mean_rt_sentence*1000, std(rt_correct(strcmp(cond_correct, 'sentence')))*1000);
                'Nonword (N)', sprintf('%6.1f ± %.1f', mean_rt_nonword*1000, std(rt_correct(strcmp(cond_correct, 'nonword')))*1000)
            };
            % Validate all entries are strings
            for i = 1:size(tableContent, 1)
                for j = 1:size(tableContent, 2)
                    if ~ischar(tableContent{i,j}) && ~isstring(tableContent{i,j})
                        tableContent{i,j} = char(string(tableContent{i,j}));
                    end
                end
            end
            tbl = Table(tableContent);
            tbl.Style = {Border('solid'), ColSep('solid'), RowSep('solid')};
            tbl.TableEntriesStyle = {HAlign('center')};
            add(rpt, tbl);
        
            % Reaction Time Distribution Plot
            add(rpt, Heading2('Reaction Time Distributions'));
            debugMode = false;
            if debugMode
                f = figure('Visible', 'on', 'Position', [100 100 800 600]);
            else
                f = figure('Visible', 'off', 'Position', [100 100 800 600]);
            end
            subplot(1,2,1)
            histogram(rt_correct(strcmp(cond_correct, 'sentence'))*1000, 'BinWidth', 50)
            title('Sentence Condition RT Distribution')
            xlabel('Reaction Time (ms)')
            ylabel('Frequency')
        
            subplot(1,2,2)
            histogram(rt_correct(strcmp(cond_correct, 'nonword'))*1000, 'BinWidth', 50)
            title('Nonword Condition RT Distribution')
            xlabel('Reaction Time (ms)')
            ylabel('Frequency')
        
            % Add the figure to the PDF
            imgObj = addImageToReport(tempFolder, f);
            add(rpt, imgObj);
            add(rpt, PageBreak());

            
            % High Gamma Plot Generation
            add(rpt, Heading2('High Gamma Plots'));
            conditionImages = high_gamma_plot_langloc(obj);

            for i = 1:size(conditionImages, 1)
                settingName = conditionImages{i, 1};
                add(rpt, Chapter(settingName));
                imageFiles = conditionImages{i,3};

                for imgPathId = 1:length(imageFiles)
                    imgPath = imageFiles{imgPathId};
                    imgObj = Image(imgPath);
                    imgObj.Width = "6in";
                    imgObj.Height = "7in";
                    add(rpt, imgObj);
                    add(rpt, PageBreak());
                end
                deleteTemporaryFiles(imageFiles)
            end
            
        case {'LangLocAudio','LangLocAudio-2'}
            mean_rt_all = mean(rt_correct, 'omitnan');
            mean_rt_sentence = mean(rt_correct(strcmp(cond_correct, 'sentence')), 'omitnan');
            mean_rt_nonword = mean(rt_correct(strcmp(cond_correct, 'nonword')), 'omitnan');
        
            % Results Table
            tableContent = {
                'Condition', 'Mean RT (ms)';
                'Overall', sprintf('%6.1f ± %.1f', mean_rt_all*1000, std(rt_correct)*1000);
                'Sentence (S)', sprintf('%6.1f ± %.1f', mean_rt_sentence*1000, std(rt_correct(strcmp(cond_correct, 'sentence')))*1000);
                'Nonword (N)', sprintf('%6.1f ± %.1f', mean_rt_nonword*1000, std(rt_correct(strcmp(cond_correct, 'nonword')))*1000)
            };
            tbl = Table(tableContent);
            tbl.Style = {Border('solid'), ColSep('solid'), RowSep('solid')};
            tbl.TableEntriesStyle = {HAlign('center')};
            add(rpt, tbl);
        
            % Reaction Time Distribution Plot
            add(rpt, Heading2('Reaction Time Distributions'));
            debugMode = false;
            if debugMode
                f = figure('Visible', 'on', 'Position', [100 100 800 600]);
            else
                f = figure('Visible', 'off', 'Position', [100 100 800 600]);
            end
            subplot(1,2,1)
            histogram(rt_correct(strcmp(cond_correct, 'sentence'))*1000, 'BinWidth', 50)
            title('Sentence Condition RT Distribution')
            xlabel('Reaction Time (ms)')
            ylabel('Frequency')
        
            subplot(1,2,2)
            histogram(rt_correct(strcmp(cond_correct, 'nonword'))*1000, 'BinWidth', 50)
            title('Nonword Condition RT Distribution')
            xlabel('Reaction Time (ms)')
            ylabel('Frequency')
        
            % Add the figure to the PDF
            imgObj = addImageToReport(tempFolder, f);
            add(rpt, imgObj);
            add(rpt, PageBreak());
        
            % Audio Duration Histogram
            add(rpt, Heading2('Audio Duration Distribution'));
            f = figure('Visible', 'off', 'Position', [100 100 800 600]);
            audioDur = obj.events_table.audio_ended_natus - obj.events_table.audio_onset_natus;
            sentence_durations = audioDur(strcmp(obj.condition, 'sentence'));
            nonword_durations = audioDur(strcmp(obj.condition, 'nonword'));
            
            subplot(1,2,1)
            histogram(sentence_durations, 'BinWidth', 0.1)
            title('Sentence Audio Duration Distribution')
            xlabel('Duration (s)')
            ylabel('Frequency')
        
            subplot(1,2,2)
            histogram(nonword_durations, 'BinWidth', 0.1)
            title('Nonword Audio Duration Distribution')
            xlabel('Duration (s)')
            ylabel('Frequency')
        
            % Add the figure to the PDF
            imgObj = addImageToReport(tempFolder, f);
            add(rpt, imgObj);
            add(rpt, PageBreak());
             tempFiles = {};
            % ================================================================
            % NEW: LangLoc responsive Electrodes – Dysoc-style visualization
            % ================================================================
            add(rpt, Chapter('LangLoc Responsive Electrodes'));
            add(rpt, Paragraph(['Electrodes with significant language responsiveness ' ...
                '(Sentences > Nonwords, cluster-corrected p < 0.05). ' ...
                'Layout: split-half scatter (top-left), mean HG per condition (top-right), ' ...
                'full time series (bottom).']));
            
            [selImgs_uni, selImgs_bip] = plot_langloc_responsive_dysoc(obj, tempFolder);
            
            if ~isempty(selImgs_uni)
                add(rpt, Heading2('Unipolar Responsive Electrodes'));
                for k = 1:numel(selImgs_uni)
                    imgObj = Image(selImgs_uni{k});
                    imgObj.Width  = '6in';
                    imgObj.Height = '7.5in';
                    add(rpt, imgObj);
                    add(rpt, PageBreak());
                end
               % deleteTemporaryFiles(selImgs_uni);
                 tempFiles = [tempFiles, selImgs_uni];
            else
                add(rpt, Paragraph('No responsive unipolar electrodes found.'));
            end
            
            if ~isempty(selImgs_bip)
                add(rpt, Heading2('Bipolar Responsive Electrodes'));
                for k = 1:numel(selImgs_bip)
                    imgObj = Image(selImgs_bip{k});
                    imgObj.Width  = '6in';
                    imgObj.Height = '7.5in';
                    add(rpt, imgObj);
                    add(rpt, PageBreak());
                end
               % deleteTemporaryFiles(selImgs_bip);
               tempFiles = [tempFiles, selImgs_bip];
            end
            % ================================================================


        
            % High Gamma Plot Generation
            add(rpt, Heading2('High Gamma Plots'));
            conditionImages = high_gamma_plot_word_boundaries_langloc(obj);

           

            for i = 1:size(conditionImages, 1)
                settingName = conditionImages{i, 1};
                add(rpt, Chapter(settingName));
                imageFiles = conditionImages{i,3};

                for imgPathId = 1:length(imageFiles)
                    imgPath = imageFiles{imgPathId};
                    imgObj = Image(imgPath);
                    imgObj.Width = "6in";
                    imgObj.Height = "7in";
                    add(rpt, imgObj);
                    add(rpt, PageBreak());
                end
                %deleteTemporaryFiles(imageFiles)
                 tempFiles = [tempFiles, imageFiles];
            end
    end

    % % Add summary report of significant channels
    % add(rpt, Chapter('Summary of Significant Channels'));
    % 
    % % Unipolar channels
    % add(rpt, Heading2('Unipolar Channels with Significant Time Clusters'));
    % sigUnipolarChannels = find(cellfun(@(x) any(x.h_sig_05), obj.stats.time_series.pSigChan));
    % if ~isempty(sigUnipolarChannels)
    %     unipolarList = cell(length(sigUnipolarChannels), 1);
    %     for i = 1:length(sigUnipolarChannels)
    %         unipolarList{i} = obj.elec_ch_label{sigUnipolarChannels(i)};
    %     end
    %     add(rpt, UnorderedList(unipolarList));
    %     add(rpt, PageBreak());
    % else
    %     add(rpt, Paragraph('No significant unipolar channels found.'));
    % end
    % 
    % % Bipolar channels
    % if(isfield(obj.stats.time_series,'pSigChan_bip'))
    %     add(rpt, Heading2('Bipolar Channels with Significant Time Clusters'));
    %     sigBipolarChannels = find(cellfun(@(x) any(x.h_sig_05), obj.stats.time_series.pSigChan_bip));
    %     if ~isempty(sigBipolarChannels)
    %         bipolarList = cell(length(sigBipolarChannels), 1);
    %         for i = 1:length(sigBipolarChannels)
    %             bipolarList{i} = obj.bip_ch_label{sigBipolarChannels(i)};
    %         end
    %         add(rpt, UnorderedList(bipolarList));
    %         add(rpt, PageBreak());
    %     else
    %         add(rpt, Paragraph('No significant bipolar channels found.'));
    %     end
    % end


     % Add summary report of significant channels
    add(rpt, Chapter('Summary of Significant LangLoc Channels'));

    % Unipolar channels
    add(rpt, Heading2('Unipolar Channels with Significant Time Clusters'));
    sigChannels = sum(cellfun(@(x) any(cell2mat(arrayfun(@(y) any(y.h_sig_05), x, 'UniformOutput', false))), (obj.stats.time_series.pSigChan_wordboundaries_langloc)'));
    sigUnipolarChannels = find(sigChannels>0);
    if ~isempty(sigUnipolarChannels)
        unipolarList = cell(length(sigUnipolarChannels), 1);
        for i = 1:length(sigUnipolarChannels)
            unipolarList{i} = obj.elec_ch_label{sigUnipolarChannels(i)};
        end
        add(rpt, UnorderedList(unipolarList));
        add(rpt, PageBreak());
    else
        add(rpt, Paragraph('No significant langloc unipolar channels found.'));
    end

    % Bipolar channels
    if(isfield(obj.stats.time_series,'pSigChan_bip'))
        add(rpt, Heading2('Bipolar Channels with Significant Time Clusters'));
        sigChannels = sum(cellfun(@(x) any(cell2mat(arrayfun(@(y) any(y.h_sig_05), x, 'UniformOutput', false))), (obj.stats.time_series.pSigChan_bip_wordboundaries_langloc)'));
        sigBipolarChannels = find(sigChannels>0);
        if ~isempty(sigBipolarChannels)
            bipolarList = cell(length(sigBipolarChannels), 1);
            for i = 1:length(sigBipolarChannels)
                bipolarList{i} = obj.bip_ch_label{sigBipolarChannels(i)};
            end
            add(rpt, UnorderedList(bipolarList));
            add(rpt, PageBreak());
        else
            add(rpt, Paragraph('No significant langloc bipolar channels found.'));
        end
    end

    % Close the PDF document   
    try
        
        close(rpt);
        deleteTemporaryFiles(tempFiles);
    catch ME
        fprintf('Error closing the report: %s\n', ME.message);
        try
            close(rpt);
            deleteTemporaryFiles(tempFiles);
        catch ME2
             deleteTemporaryFiles(tempFiles);
            fprintf('Failed to close the report again: %s\n', ME2.message);
        end
    end

    % Save the updated object in the crunched folder
    saveUpdatedObject(obj);
end

% -----------------------------------------------------------------------
function [imageFiles_uni, imageFiles_bip] = plot_langloc_responsive_dysoc(obj, tempFolder)
% Dysoc-style per-electrode plots for langloc-responsive channels (S > N).
% Uses word-boundary concatenation to handle variable stimulus durations.

epochTimeRange = [-0.5 0.5];
numWords       = 12;
colors         = [0 0.4470 0.7410; 0.8500 0.3250 0.0980];  % blue / red
conds          = obj.condition;

if contains(obj.experiment, 'MITLangloc')
    conditions = {'Sentences', 'Jabberwocky'};
    condLabels = {'Sentences', 'Jabberwocky'};
else
    conditions = {'sentence', 'nonword'};
    condLabels = {'Sentences', 'Nonwords'};
end

% ---- Extract per-word epochs and concatenate (mirrors word_boundaries function) ----
fprintf('plot_langloc_responsive_dysoc: extracting per-word epochs...\n');
concatS = []; concatN = []; concatS_bip = []; concatN_bip = [];
epochData_bip_ref = [];

for wordPos = 1:numWords
    fprintf('  Word position %d/%d...\n', wordPos, numWords);
    [epochData, epochData_bip] = obj.extract_trial_epochs( ...
        'epoch_tw', epochTimeRange, 'probe_key', wordPos + 1);

    concatS = cat(3, concatS, epochData(:, ismember(conds, conditions{1}), :));
    concatN = cat(3, concatN, epochData(:, ismember(conds, conditions{2}), :));

    if ~isempty(epochData_bip)
        concatS_bip = cat(3, concatS_bip, epochData_bip(:, ismember(conds, conditions{1}), :));
        concatN_bip = cat(3, concatN_bip, epochData_bip(:, ismember(conds, conditions{2}), :));
        epochData_bip_ref = epochData_bip;
    end
end

% Word-boundary geometry
timePointsPerWord = size(concatS, 3) / numWords;
totalTimePoints   = timePointsPerWord * numWords;
wordBoundaries    = 0 : timePointsPerWord : totalTimePoints;
x_plot            = 1:totalTimePoints;
midPoints         = (wordBoundaries(1:end-1) + wordBoundaries(2:end)) / 2;
wordLabels        = arrayfun(@(w) sprintf('W%d', w), 1:numWords, 'UniformOutput', false);

% ---- Compute S>N stats on concatenated data: unipolar ----
fprintf('  Computing S>N selectivity stats (unipolar)...\n');
nChan_uni = size(concatS, 1);
pSig_uni  = cell(1, nChan_uni);
parfor iChan = 1:nChan_uni
    pSig_uni{iChan} = timePermCluster( ...
        squeeze(concatS(iChan,:,:)), squeeze(concatN(iChan,:,:)), 'numTail', 1);
end
responsiveIdx_uni = find(cellfun(@(s) any(s.h_sig_05), pSig_uni));
fprintf('  Found %d responsive unipolar electrodes.\n', numel(responsiveIdx_uni));

% ---- Compute S>N stats: bipolar ----
pSig_bip = {}; responsiveIdx_bip = []; chanLab_bip = {};
if ~isempty(concatS_bip)
    chanLab_bip = obj.bip_ch_label;
    nChan_bip   = size(concatS_bip, 1);
    pSig_bip    = cell(1, nChan_bip);
    parfor iChan = 1:nChan_bip
        pSig_bip{iChan} = timePermCluster( ...
            squeeze(concatS_bip(iChan,:,:)), squeeze(concatN_bip(iChan,:,:)), 'numTail', 1);
    end
    responsiveIdx_bip = find(cellfun(@(s) any(s.h_sig_05), pSig_bip));
    fprintf('  Found %d responsive bipolar electrodes.\n', numel(responsiveIdx_bip));
end

chanLab_uni = obj.elec_ch_label(obj.elec_ch_valid);

% Store stats for optional reuse downstream
obj.stats.time_series.pSigChan_langloc_responsive_uni = pSig_uni;
if ~isempty(pSig_bip)
    obj.stats.time_series.pSigChan_langloc_responsive_bip = pSig_bip;
end

% Pack plot geometry into a struct for clean argument passing
P.wordBoundaries    = wordBoundaries;
P.timePointsPerWord = timePointsPerWord;
P.numWords          = numWords;
P.totalTimePoints   = totalTimePoints;
P.x_plot            = x_plot;
P.midPoints         = midPoints;
P.wordLabels        = wordLabels;
P.colors            = colors;
P.condLabels        = condLabels;

% ---- Generate figures ----
imageFiles_uni = generate_responsive_figs_wordboundary( ...
    concatS, concatN, pSig_uni, responsiveIdx_uni, chanLab_uni, P, 'unipolar', tempFolder);

imageFiles_bip = {};
if ~isempty(concatS_bip)
    imageFiles_bip = generate_responsive_figs_wordboundary( ...
        concatS_bip, concatN_bip, pSig_bip, responsiveIdx_bip, chanLab_bip, P, 'bipolar', tempFolder);
end
end

% -----------------------------------------------------------------------
function imageFiles = generate_responsive_figs_wordboundary( ...
    dataSentence, dataNonword, pSig, responsiveIdx, chanLab, P, dataType, tempFolder)
% Renders one dysoc-style figure per responsive electrode using word-boundary layout.

imageFiles = {};

for k = 1:numel(responsiveIdx)
    iChan     = responsiveIdx(k);
    chanLabel = chanLab{iChan};

    % Safely squeeze channel data -> [nTrials x nTime]
    dataS = reshape(dataSentence(iChan,:,:), size(dataSentence,2), size(dataSentence,3));
    dataN = reshape(dataNonword(iChan,:,:),  size(dataNonword,2),  size(dataNonword,3));

    f = figure('Visible', 'off', 'Position', [100 100 1400 900], 'Renderer', 'painters');
    sgtitle(sprintf('electrode id: %d  (%s)', iChan, chanLabel), ...
        'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'none');

    % ---- TOP-LEFT: split-half scatter ----
    ax1  = subplot(2, 2, 1);
    allX = []; allY = []; allC = [];

    for iCond = 1:2
        D = [dataS; dataN * 0];  % placeholder
        if iCond == 1, D = dataS; else, D = dataN; end
        n = size(D, 1);
        if n < 2, continue; end

        rng('default');
        perm  = randperm(n);
        h1    = perm(1 : floor(n/2));
        h2    = perm(floor(n/2)+1 : 2*floor(n/2));
        nPair = numel(h1);

        % Per-trial mean HG across all word positions
        hg1 = mean(D(h1, :), 2);
        hg2 = mean(D(h2, :), 2);
        allX = [allX; hg1];
        allY = [allY; hg2];
        allC = [allC; repmat(P.colors(iCond,:), nPair, 1)];
    end

    if numel(allX) > 1
        [r_val, p_val] = corr(allX, allY);
        n_total        = numel(allX) * 2;
    else
        r_val = NaN; p_val = NaN; n_total = 0;
    end

    scatter(ax1, allX, allY, 30, allC, 'filled', 'MarkerFaceAlpha', 0.7);
    hold(ax1, 'on');
    allVals = [allX; allY];
    lo = min(allVals) - 0.05;  hi = max(allVals) + 0.05;
    plot(ax1, [lo hi], [lo hi], 'k-', 'LineWidth', 1);
    xlim(ax1, [lo hi]);  ylim(ax1, [lo hi]);
    title(ax1, sprintf('splithalf r = %.2f, p = %.2f\n(n = %d), avg of 1 splits', ...
        r_val, p_val, n_total), 'FontSize', 9);
    xlabel(ax1, 'HG power');  ylabel(ax1, 'HG power');
    hold(ax1, 'off');

    % ---- TOP-RIGHT: bar chart – one bar per condition, mean across all word positions ----
    ax2 = subplot(2, 2, 2);

    meanPerTrialS = mean(dataS, 2);   % [nTrialsS x 1]
    meanPerTrialN = mean(dataN, 2);   % [nTrialsN x 1]
    condMeans     = [mean(meanPerTrialS), mean(meanPerTrialN)];
    condSEMs      = [std(meanPerTrialS)/sqrt(size(dataS,1)), ...
                     std(meanPerTrialN)/sqrt(size(dataN,1))];

    b = bar(ax2, 1:2, condMeans, 'FaceColor', 'flat');
    b.CData(1,:) = P.colors(1,:);
    b.CData(2,:) = P.colors(2,:);
    hold(ax2, 'on');
    errorbar(ax2, 1:2, condMeans, condSEMs, 'k', ...
        'LineStyle', 'none', 'LineWidth', 1.5, 'CapSize', 8);
    set(ax2, 'XTick', 1:2, 'XTickLabel', P.condLabels, 'FontSize', 10);
    ylabel(ax2, 'HG power');
    title(ax2, 'Mean HG Power (all word positions)');
    barTop = max(condMeans + condSEMs);
    if barTop > 0, ylim(ax2, [0, barTop * 1.35]); end
    hold(ax2, 'off');

    % ---- BOTTOM: word-boundary time series (spans both columns) ----
    ax3 = subplot(2, 1, 2);
    hold(ax3, 'on');

    mnS = nanmean(dataS, 1);
    seS = nanstd(dataS, 0, 1) / sqrt(size(dataS, 1));
    mnN = nanmean(dataN, 1);
    seN = nanstd(dataN, 0, 1) / sqrt(size(dataN, 1));

    patch(ax3, [P.x_plot fliplr(P.x_plot)], [mnS+seS fliplr(mnS-seS)], ...
        P.colors(1,:), 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    hS = plot(ax3, P.x_plot, mnS, 'Color', P.colors(1,:), 'LineWidth', 2);

    patch(ax3, [P.x_plot fliplr(P.x_plot)], [mnN+seN fliplr(mnN-seN)], ...
        P.colors(2,:), 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    hN = plot(ax3, P.x_plot, mnN, 'Color', P.colors(2,:), 'LineWidth', 2);

    % Word boundary markers
    for bnd = P.wordBoundaries
        if bnd > 0 && bnd < P.totalTimePoints
            xline(ax3, bnd, 'k--', 'LineWidth', 1.5, 'Alpha', 0.6);
        end
    end

    % Significance dots (S > N, red) below traces
    sigT = find(pSig{iChan}.h_sig_05 == 1);
    yL   = ylim(ax3);
    if ~isempty(sigT)
        scatter(ax3, sigT, (yL(1) - 0.05) * ones(size(sigT)), 10, 'r', 'filled');
    end

    % Word-position x-axis labels (W1, W2, ...)
    set(ax3, 'XTick', P.midPoints, 'XTickLabel', P.wordLabels, 'FontSize', 9);
    yline(ax3, 0, 'k-', 'LineWidth', 0.5, 'Alpha', 0.4);
    xlabel(ax3, 'Word Position');
    ylabel(ax3, 'HG power (gamma filtered signal)');
    legend(ax3, [hS hN], P.condLabels, 'Location', 'best');
    xlim(ax3, [P.x_plot(1) P.x_plot(end)]);
    hold(ax3, 'off');

    % Save
    safeName = regexprep(chanLabel, '[^a-zA-Z0-9]', '_');
    fname    = fullfile(tempFolder, sprintf('responsive_%s_%d_%s.png', dataType, iChan, safeName));
    exportgraphics(f, fname, 'Resolution', 300);
    imageFiles{end+1} = fname;
    close(f);
end
end




function conditionImages = high_gamma_plot_langloc(obj)
    import mlreportgen.report.*
    import mlreportgen.dom.*

    % Temporary folder for saving images
    tempFolder = tempdir;

    acc = [obj.events_table.accuracy];
    numTrials = length(obj.condition);
    conditionIds = obj.condition;

    % Define settings specific to langloc experiments
    conditionImages = {
        'All Trials', 1:numTrials;
    };
    
    accurateTrials = find(acc == 1);
    if ~isempty(accurateTrials)
        conditionImages(end+1, :) = {'Accurate Trials', accurateTrials};
    end

    % Data epoching - use appropriate time range for langloc
    epochTimeRange = [-0.5 6];
    [epochData, epochData_bip] = obj.extract_trial_epochs('epoch_tw', epochTimeRange, 'selectChannels', obj.elec_ch_clean);

    for i = 1:size(conditionImages, 1)
        settingName = conditionImages{i, 1};
        trials2include = conditionImages{i, 2};
        imageFiles = {};

        % Process unipolar data
        unipolarImages = process_and_save_images_langloc(obj, epochData, conditionIds, trials2include, 'unipolar', obj.elec_ch_label(obj.elec_ch_valid), epochTimeRange, tempFolder, settingName);
        imageFiles = [imageFiles unipolarImages];

        % Process bipolar data
        if ~isempty(epochData_bip)
            bipolarImages = process_and_save_images_langloc(obj, epochData_bip, conditionIds, trials2include, 'bipolar', obj.bip_ch_label, epochTimeRange, tempFolder, settingName);
            imageFiles = [imageFiles bipolarImages];
        end

        conditionImages{i,3} = imageFiles;
    end
end
function imageFiles = process_and_save_images_langloc(obj, data, conds, trials2include, data_type, chanLab, epochTimeRange, tempFolder, settingName)
    duration = size(data, 3);
    x = linspace(epochTimeRange(1), epochTimeRange(2), duration);

    data2process = data(:,trials2include,:);
    conds2process = conds(trials2include);
    
    numChan = 10; % Number of channels to plot per figure
    totChanBlock = ceil(size(data2process, 1) / numChan);

    % Define langloc-specific conditions
    if contains(obj.experiment, {'LangLoc', 'LangLocVisual'}) || contains(obj.experiment, 'MITLangloc')
        if contains(obj.experiment, 'MITLangloc')
            conditions = {'Sentences', 'Jabberwocky'};
        else
            conditions = {'sentence', 'nonword'};
        end
    else % LangLocAudio experiments
        conditions = {'sentence', 'nonword'};
    end

    colors = [0 0.4470 0.7410; 0.8500 0.3250 0.0980]; % Blue and red for sentence/nonword

    imageFiles = {};
    
    assert(size(data2process,2)==length(conds2process),'Uneven trial length');
    
    % Statistical analysis for two conditions (sentence vs nonword)
    aCondition = conditions{1};
    bCondition = conditions{2};

    if (~ismember(aCondition, conds2process)) || (~ismember(bCondition, conds2process))
        warning('Specified conditions "%s" and "%s" do not match available conditions for %s.', aCondition, bCondition, obj.experiment);
        return;
    end

    aTrialData = squeeze(data2process(:, ismember(conds2process, aCondition), :));
    bTrialData = squeeze(data2process(:, ismember(conds2process, bCondition), :));

    % Use timePermCluster for sentence vs nonword comparison
    pSig = cell(1, size(data2process, 1));
    parfor iChan = 1:size(data2process, 1)
        pSig{iChan} = timePermCluster(squeeze(aTrialData(iChan, :, :)), squeeze(bTrialData(iChan, :, :)), 'numTail', 1);
    end
    
    % Save pSig results
    if strcmp(data_type, 'unipolar')
        obj.stats.time_series.(['pSigChan_langloc_' strrep(settingName, ' ', '_')]) = pSig;
    else
        obj.stats.time_series.(['pSigChan_bip_langloc_' strrep(settingName, ' ', '_')]) = pSig;
    end

    % Generate plots
    for iF = 0:totChanBlock-1
        f = figure('Visible', 'off', 'Position', [100 100 1000 1200], 'Renderer', 'painters');
        
        numRows = 5;
        numCols = 2;

        for iChan = 1:min(numChan, size(data2process, 1) - iF*numChan)
            iChan2 = iF*numChan + iChan;

            subplot(numRows, numCols, iChan);
            hold on;

            % Get significant time points from original HG analysis
            if strcmp(data_type, 'unipolar')
                sigTime = x(obj.stats.time_series.pSigChan{iChan2}.h_sig_05 == 1);
            else
                sigTime = x(obj.stats.time_series.pSigChan_bip{iChan2}.h_sig_05 == 1);
            end

            % Plot data for each condition
            legendEntries = {};
            for iCond = 1:length(conditions)
                condTrials = ismember(conds2process, conditions{iCond});
                condData = squeeze(data2process(iChan2, condTrials, :));
                
                condMean = nanmean(condData, 1);
                condSEM = nanstd(condData, 0, 1) / sqrt(size(condData, 1));
                
                % Plot with condition-specific color
                patch([x fliplr(x)], [condMean+condSEM fliplr(condMean-condSEM)], ...
                      colors(iCond,:), 'FaceAlpha', 0.3, 'EdgeColor', 'none');
                plot(x, condMean, 'Color', colors(iCond,:), 'LineWidth', 1.5);
                
                legendEntries{end+1} = conditions{iCond};
            end

            % Calculate appropriate y-limits
            allData = squeeze(data2process(iChan2, :, :));
            allMeans = nanmean(allData, 1);
            allSEMs = nanstd(allData, 0, 1) / sqrt(size(allData, 1));
            max_val = max(allMeans + allSEMs);
            ylim([-1 max_val + 0.5]);

            % Add significance markers
            scatter(sigTime, -0.5 * ones(size(sigTime)), 10, 'k', 'filled');
            
            % Add condition contrast significance
            sigTimeCond = x(pSig{iChan2}.h_sig_05 == 1);
            scatter(sigTimeCond, -0.75 * ones(size(sigTimeCond)), 10, 'r', 'filled');

            % Add baseline and formatting
            plot([epochTimeRange(1) epochTimeRange(2)], [0 0], 'k', 'LineWidth', 0.25);
            plot([0 0], ylim, 'k--', 'LineWidth', 0.5); % Stimulus onset
            
            title(chanLab{iChan2}, 'Interpreter', 'none', 'FontSize', 12);
            xlabel('Time (s)', 'FontSize', 10);
            ylabel('HG Power', 'FontSize', 6);
            set(gca, 'FontSize', 12);
            
            hold off;
        end

        % Add comprehensive legend
        legendAxes = axes('Position', [0.1, 0.02, 0.8, 0.05], 'Visible', 'off');
        hold(legendAxes, 'on');
        
        % Legend elements
        scatter(legendAxes, NaN, NaN, 30, 'k', 'filled');
        scatter(legendAxes, NaN, NaN, 30, 'r', 'filled');
        for iCond = 1:length(conditions)
            plot(legendAxes, NaN, NaN, 'Color', colors(iCond,:), 'LineWidth', 2);
        end
        
        legendLabels = {'Significant HG activations', 'Significant Sentence vs Nonword'};
        for iCond = 1:length(conditions)
            legendLabels{end+1} = conditions{iCond};
        end
        
        legend(legendAxes, legendLabels, 'Orientation', 'horizontal', 'Location', 'southoutside');
        hold(legendAxes, 'off');

        % Save figure
        imageFileName = fullfile(tempFolder, sprintf('%s_langloc_%s_Block%d.png', data_type, strrep(settingName, ' ', '_'), iF));
        exportgraphics(f, imageFileName, 'Resolution', 300);
        imageFiles{end+1} = imageFileName;
        
        close(f);
    end
end

function conditionImages = high_gamma_plot_word_boundaries_langloc(obj)
    import mlreportgen.report.*;
    import mlreportgen.dom.*;

    % Verify this is a langloc experiment
    langlocExperiments = {'LangLocVisual', 'LangLoc', 'MITLangloc', 'LangLocAudio', 'LangLocAudio-2'};
    if ~ismember(obj.experiment, langlocExperiments)
        warning('Word boundaries analysis is designed for langloc experiments. Current experiment: %s', obj.experiment);
        conditionImages = {};
        return;
    end

    % Temporary folder for saving images
    tempFolder = tempdir;
    if ~exist(tempFolder, 'dir')
        mkdir(tempFolder);
    end

    % Define settings specific to langloc experiments
    acc = [obj.events_table.accuracy];
    numTrials = size(obj.trial_timing, 1);
    
    conditionImages = {
        'All Trials', 1:numTrials;
    };
    
    % Add accurate trials condition
    accurateTrials = find(acc == 1);
    if ~isempty(accurateTrials)
        conditionImages(end+1, :) = {'Accurate Trials', accurateTrials};
    end

    % Data epoching parameters
    epochTimeRange = [-0.5 0.5];
    numWords = 12;
    
    % Validate that we have the required conditions
    availableConditions = unique(obj.condition);
    if ~ismember('sentence', availableConditions) || ~ismember('nonword', availableConditions)
        error('Required conditions "sentence" and "nonword" not found. Available conditions: %s', strjoin(availableConditions, ', '));
    end

    % Initialize statistical analysis storage
    pSig = cell(size(obj.elec_ch_label, 1), numWords);
    pSig_bip = cell(size(obj.bip_ch_label, 1), numWords);
    
    % Loop through each word position to extract and analyze data
    concatenatedEpochsSentence = [];
    concatenatedEpochsNonword = [];
    concatenatedEpochsSentence_bip = [];
    concatenatedEpochsNonword_bip = [];
    
    fprintf('Processing word boundaries for %s experiment...\n', obj.experiment);
    
    for wordPos = 1:numWords
        fprintf('  Processing word position %d/%d...\n', wordPos, numWords);
        
        % Extract epochs for current word position
        [epochData, epochData_bip] = obj.extract_trial_epochs('epoch_tw', epochTimeRange, 'probe_key', wordPos+1);

        % Extract sentence and nonword trials
        sentenceTrials = epochData(:, ismember(obj.condition, 'sentence'), :);
        nonwordTrials = epochData(:, ismember(obj.condition, 'nonword'), :);

        % Concatenate along time axis
        if isempty(concatenatedEpochsSentence)
            concatenatedEpochsSentence = sentenceTrials;
            concatenatedEpochsNonword = nonwordTrials;
        else
            concatenatedEpochsSentence = cat(3, concatenatedEpochsSentence, sentenceTrials);
            concatenatedEpochsNonword = cat(3, concatenatedEpochsNonword, nonwordTrials);
        end

        % Perform timePermCluster test for significance (unipolar)
        fprintf('Performing cluster stat ..')
        parfor iChan = 1:size(epochData, 1)

            aTrialData = squeeze(epochData(iChan, ismember(obj.condition, 'sentence'), :));
            bTrialData = squeeze(epochData(iChan, ismember(obj.condition, 'nonword'), :));
            
            if size(aTrialData, 1) > 1 && size(bTrialData, 1) > 1
                pSig{iChan, wordPos} = timePermCluster(aTrialData, bTrialData, 'numTail', 1);
            else
                % Handle insufficient trials
                pSig{iChan, wordPos} = struct('h_sig_05', zeros(1, size(epochData, 3)), 'p_val', ones(1, size(epochData, 3)));
            end
            fprintf(1,'.');
        end
        fprintf(1,'] done\n')
        % Perform timePermCluster test for bipolar data
        if ~isempty(epochData_bip)
            % Handle bipolar concatenation
            sentenceTrials_bip = epochData_bip(:, ismember(obj.condition, 'sentence'), :);
            nonwordTrials_bip = epochData_bip(:, ismember(obj.condition, 'nonword'), :);
            
            if isempty(concatenatedEpochsSentence_bip)
                concatenatedEpochsSentence_bip = sentenceTrials_bip;
                concatenatedEpochsNonword_bip = nonwordTrials_bip;
            else
                concatenatedEpochsSentence_bip = cat(3, concatenatedEpochsSentence_bip, sentenceTrials_bip);
                concatenatedEpochsNonword_bip = cat(3, concatenatedEpochsNonword_bip, nonwordTrials_bip);
            end
            
            parfor iChan = 1:size(epochData_bip, 1)
                aTrialData = squeeze(epochData_bip(iChan, ismember(obj.condition, 'sentence'), :));
                bTrialData = squeeze(epochData_bip(iChan, ismember(obj.condition, 'nonword'), :));
                
                if size(aTrialData, 1) > 1 && size(bTrialData, 1) > 1
                    pSig_bip{iChan, wordPos} = timePermCluster(aTrialData, bTrialData, 'numTail', 1);
                else
                    % Handle insufficient trials
                    pSig_bip{iChan, wordPos} = struct('h_sig_05', zeros(1, size(epochData_bip, 3)), 'p_val', ones(1, size(epochData_bip, 3)));
                end
            end
        end
    end

    % Calculate word boundary parameters
    timePointsPerWord = size(concatenatedEpochsSentence, 3) / numWords;
    totalTimePoints = timePointsPerWord * numWords;
    wordBoundaries = 0:timePointsPerWord:totalTimePoints;

    % Store statistical results
    obj.stats.time_series.pSigChan_wordboundaries_langloc = pSig;
    if ~isempty(epochData_bip)
        obj.stats.time_series.pSigChan_bip_wordboundaries_langloc = pSig_bip;
    end

    % Process and plot data for each condition setting
    for i = 1:size(conditionImages, 1)
        settingName = conditionImages{i, 1};
        trials2include = conditionImages{i, 2};
        imageFiles = {};

        fprintf('  Generating plots for %s...\n', settingName);

        % Process unipolar data
        imageFiles_unipolar = process_and_save_images_word_boundaries_langloc(obj, ...
            concatenatedEpochsSentence, concatenatedEpochsNonword, pSig, wordBoundaries, ...
            timePointsPerWord, totalTimePoints, tempFolder, settingName, 'unipolar', obj.elec_ch_label);
        imageFiles = [imageFiles imageFiles_unipolar];

        % Process bipolar data if available
        if ~isempty(epochData_bip)
            imageFiles_bipolar = process_and_save_images_word_boundaries_langloc(obj, ...
                concatenatedEpochsSentence_bip, concatenatedEpochsNonword_bip, pSig_bip, wordBoundaries, ...
                timePointsPerWord, totalTimePoints, tempFolder, settingName, 'bipolar', obj.bip_ch_label);
            imageFiles = [imageFiles imageFiles_bipolar];
        end

        % Store image files for this setting
        conditionImages{i, 3} = imageFiles;
    end
    
    fprintf('Word boundaries analysis completed for %s\n', obj.experiment);
end

function imageFiles = process_and_save_images_word_boundaries_langloc(obj, ...
    dataSentence, dataNonword, pSig, wordBoundaries, timePointsPerWord, ...
    totalTimePoints, tempFolder, settingName, dataType, chanLab)
    
    numChanBlock = 5; % Number of channels per figure
    totChanBlock = ceil(size(dataSentence, 1) / numChanBlock);
    
    % Langloc-specific color scheme
    colors = [0 0.4470 0.7410; 0.8500 0.3250 0.0980]; % Blue for sentences, red for nonwords
    
    imageFiles = {};
    
    % Calculate number of words for processing
    numWords = size(dataSentence, 3) / timePointsPerWord;
    
    for iChanBlock = 0:totChanBlock-1
        f = figure('Visible', 'off', 'Position', [100 100 1200 1400], 'Renderer', 'painters');
        
        for iChan = 1:min(numChanBlock, size(dataSentence, 1) - iChanBlock*numChanBlock)
            iChan2 = iChanBlock*numChanBlock + iChan;
            
            subplot(numChanBlock, 1, iChan);
            hold on;
            
            % Enhanced title with channel information
            title(sprintf('%s ', chanLab{iChan2}), 'Interpreter', 'none', 'FontSize', 20, 'FontWeight', 'bold');
            
            % Process sentence data
            trialData_sentence = squeeze(dataSentence(iChan2, :, :));
            trialMean_sentence = nanmean(trialData_sentence, 1);
            trialSEM_sentence = nanstd(trialData_sentence, 0, 1) / sqrt(size(trialData_sentence, 1));
            max_val_sentence = max(trialMean_sentence + trialSEM_sentence);
            
            % Process nonword data
            trialData_nonword = squeeze(dataNonword(iChan2, :, :));
            trialMean_nonword = nanmean(trialData_nonword, 1);
            trialSEM_nonword = nanstd(trialData_nonword, 0, 1) / sqrt(size(trialData_nonword, 1));
            max_val_nonword = max(trialMean_nonword + trialSEM_nonword);
            
            x = 1:totalTimePoints; % Time points for plotting
            
            % Plot sentence data with enhanced visualization
             plot(x, trialMean_sentence, 'Color', colors(1,:), 'LineWidth', 2);
            patch([x fliplr(x)], [trialMean_sentence+trialSEM_sentence fliplr(trialMean_sentence-trialSEM_sentence)], ...
                  colors(1,:), 'FaceAlpha', 0.3, 'EdgeColor', 'none');
           
            
            % Plot nonword data with enhanced visualization
            plot(x, trialMean_nonword, 'Color', colors(2,:), 'LineWidth', 2);
            patch([x fliplr(x)], [trialMean_nonword+trialSEM_nonword fliplr(trialMean_nonword-trialSEM_nonword)], ...
                  colors(2,:), 'FaceAlpha', 0.3, 'EdgeColor', 'none');
            
            
            % Add enhanced word boundary markers
            for boundary = wordBoundaries
                if boundary > 0 && boundary <= totalTimePoints
                    xline(boundary, 'k--', 'LineWidth', 1.5, 'Alpha', 0.7);
                end
            end
            
            % Set x-axis labels with improved formatting
            if numWords <= 12
                midPoints = (wordBoundaries(1:end-1) + wordBoundaries(2:end)) / 2;
                wordLabels = arrayfun(@(x) sprintf('W%d', x), 1:numWords, 'UniformOutput', false);
                set(gca, 'XTick', midPoints, 'XTickLabel', wordLabels);
            else
                % For more than 12 words, show every other word
                midPoints = (wordBoundaries(1:2:end-1) + wordBoundaries(2:2:end)) / 2;
                wordLabels = arrayfun(@(x) sprintf('W%d', x), 1:2:numWords, 'UniformOutput', false);
                set(gca, 'XTick', midPoints, 'XTickLabel', wordLabels);
            end
            
            % Plot significance markers with enhanced visualization
            sigTimePoints = [];
            for wordPos = 1:numWords
                if wordPos <= size(pSig, 2) && ~isempty(pSig{iChan2, wordPos})
                    wordSigPoints = find(pSig{iChan2, wordPos}.h_sig_05) + (wordPos-1)*timePointsPerWord;
                    sigTimePoints = [sigTimePoints, wordSigPoints];
                end
            end
            
            % Determine overall max_val for proper scaling
            max_val = max(max_val_sentence, max_val_nonword);
            
            % Set enhanced y-axis limits
            ylim([-1.5 max_val + 0.5]);
            
            % Add significance markers
            if ~isempty(sigTimePoints)
                scatter(sigTimePoints, -1.25*ones(size(sigTimePoints)), 15, 'r', 'filled',  'LineWidth', 0.5);
            end
            
            % Add reference lines
            yline(0, 'k-.', 'LineWidth', 1, 'Alpha', 0.5);
            
            % Enhanced axis formatting
            xlabel('Word Position', 'FontSize', 10);
            ylabel('HG power', 'FontSize', 6);
            
            grid on;
            grid minor;
            
            % Add trial count information
            % text(0.02, 0.98, sprintf('Sentence trials: %d\nNonword trials: %d', ...
            %     size(trialData_sentence, 1), size(trialData_nonword, 1)), ...
            %     'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', 9, ...
            %     'BackgroundColor', 'white', 'EdgeColor', 'black', 'FaceAlpha', 0.8);
            
            hold off;
        end
        
        % Enhanced legend with langloc-specific information
        legendAxes = axes('Position', [0.1, 0.01, 0.8, 0.08], 'Visible', 'off');
        hold(legendAxes, 'on');
        
        % Legend elements
        scatter(legendAxes, NaN, NaN, 30, 'r', 'filled');
        plot(legendAxes, NaN, NaN, 'Color', colors(1,:), 'LineWidth', 2);
        plot(legendAxes, NaN, NaN, 'Color', colors(2,:), 'LineWidth', 2);
        plot(legendAxes, NaN, NaN, 'k--', 'LineWidth', 1.5);
        plot(legendAxes, NaN, NaN, 'k-.', 'LineWidth', 1);
        
        legendLabels = {
            'Significant Sentence vs Nonword',
            'Sentences',
            'Nonwords',
            'Word Boundaries',
            'Baseline'
        };
        
        legend(legendAxes, legendLabels, 'Orientation', 'horizontal', 'Location', 'southoutside', 'FontSize', 20);
        hold(legendAxes, 'off');
        
        % Enhanced figure title
        sgtitle(sprintf('%s - %s Word Boundaries Analysis (%s)', ...
            obj.experiment, strrep(settingName, '_', ' '), upper(dataType)), 'FontSize', 20, 'FontWeight', 'bold');

        % Save the figure with descriptive filename
        imageFileName = fullfile(tempFolder, sprintf('%s_wordboundaries_langloc_%s_%s_Block%d.png', ...
            dataType, obj.experiment, strrep(settingName, ' ', '_'), iChanBlock));
        exportgraphics(f, imageFileName, 'Resolution', 300);
        imageFiles{end+1} = imageFileName;
        
        % Close the figure to free memory
        close(f);
    end
    
    % Display summary statistics
    sigChannels = sum(cellfun(@(x) any(cell2mat(arrayfun(@(y) any(y.h_sig_05), x, 'UniformOutput', false))), pSig'));
    totalSigChannels = sum(sigChannels>0);

    fprintf('    %s analysis completed: %d channels, %d significant channels\n', ...
        upper(dataType), size(dataSentence, 1), totalSigChannels);
end


function deleteTemporaryFiles(imageFiles)
    for iFile = 1:numel(imageFiles)
        delete(imageFiles{iFile}); % Delete each temporary file
    end
end


function saveUpdatedObject(obj)
    % Create the crunched folder if it doesn't exist
    crunchedFolder = fullfile(obj.crunched_file_path);
    if ~exist(crunchedFolder, 'dir')
        mkdir(crunchedFolder);
    end

    % Generate the filename
    filename = fullfile(crunchedFolder, [obj.subject '_' obj.experiment '_crunched_HG_ZScore.mat']);

    % Save the object
    save(filename, 'obj', '-v7.3');
    
    fprintf('Updated object saved as: %s\n', filename);
end

function imgObj = addImageToReport(tempFolder, f)
    import mlreportgen.report.*
    import mlreportgen.dom.*

    % Generate the image file name
    imageFileName = fullfile(tempFolder, ['imageAdd_' datestr(now, 'yyyymmdd_HHMMSS_FFF') '.png']);
    
    % Save the figure as a PNG file
    exportgraphics(f, imageFileName,'Resolution',300);

    % Verify the file was created successfully
    if ~exist(imageFileName, 'file')
        error('Failed to create image file: %s', imageFileName);
    end
    
    % Create an Image object from the saved file
    imgObj = Image(imageFileName);
    imgObj.Width = "6in";
    imgObj.Height = "7in";  
   
    
    % Close the figure to free up memory
    close(f);
end
