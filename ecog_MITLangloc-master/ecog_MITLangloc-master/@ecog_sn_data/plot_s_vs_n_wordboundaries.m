
function plot_s_vs_n_wordboundaries(obj, ops)
    % PLOT_S_VS_N_WORDBOUNDARIES - Plot results of word boundaries time-series permutation test

    % Extract stored data
    concatenatedEpochsSentence = obj.s_vs_n_wordboundaries_concatenatedEpochsSentence;
    concatenatedEpochsNonword = obj.s_vs_n_wordboundaries_concatenatedEpochsNonword;
    pSig = obj.s_vs_n_wordboundaries_pSig;
    wordBoundaries = obj.s_vs_n_wordboundaries_wordBoundaries;
    timePointsPerWord = obj.s_vs_n_wordboundaries_timePointsPerWord;
    totalTimePoints = obj.s_vs_n_wordboundaries_totalTimePoints;

    % Plot parameters
    numChanBlock = 5; % Number of channels per figure
    totChanBlock = ceil(size(concatenatedEpochsSentence, 1) / numChanBlock);

    % Color scheme
    colors = [0 0.4470 0.7410; 0.8500 0.3250 0.0980]; % Blue for sentences, red for nonwords

    % Calculate time axis
    x = 1:totalTimePoints;
    numWords = size(pSig, 2);

    % Generate plots for unipolar channels
    for iChanBlock = 0:(totChanBlock-1)
        f = figure('Visible', 'off', 'Position', [100 100 1200 1400], 'Renderer', 'painters');

        for iChan = 1:min(numChanBlock, size(concatenatedEpochsSentence, 1) - iChanBlock*numChanBlock)
            iChan2 = iChanBlock*numChanBlock + iChan;
            subplot(numChanBlock, 1, iChan);
            hold on;

            % Enhanced title with channel information
            title(sprintf('%s', obj.elec_ch_label{iChan2}), 'Interpreter', 'none', ...
                  'FontSize', 10, 'FontWeight', 'bold');

            % Process sentence data
            trialData_sentence = squeeze(concatenatedEpochsSentence(iChan2, :, :));
            trialMean_sentence = nanmean(trialData_sentence, 1);
            trialSEM_sentence = nanstd(trialData_sentence, 0, 1) / sqrt(size(trialData_sentence, 1));
            maxval_sentence = max(trialMean_sentence + trialSEM_sentence);

            % Process nonword data
            trialData_nonword = squeeze(concatenatedEpochsNonword(iChan2, :, :));
            trialMean_nonword = nanmean(trialData_nonword, 1);
            trialSEM_nonword = nanstd(trialData_nonword, 0, 1) / sqrt(size(trialData_nonword, 1));
            maxval_nonword = max(trialMean_nonword + trialSEM_nonword);

            % Plot sentence data with enhanced visualization
            plot(x, trialMean_sentence, 'Color', colors(1,:), 'LineWidth', 2);
            patch([x, fliplr(x)], [trialMean_sentence+trialSEM_sentence, ...
                  fliplr(trialMean_sentence-trialSEM_sentence)], ...
                  colors(1,:), 'FaceAlpha', 0.3, 'EdgeColor', 'none');

            % Plot nonword data with enhanced visualization
            plot(x, trialMean_nonword, 'Color', colors(2,:), 'LineWidth', 2);
            patch([x, fliplr(x)], [trialMean_nonword+trialSEM_nonword, ...
                  fliplr(trialMean_nonword-trialSEM_nonword)], ...
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
                wordLabels = arrayfun(@(x) sprintf('Wd%d', x), 1:numWords, 'UniformOutput', false);
                set(gca, 'XTick', midPoints, 'XTickLabel', wordLabels);
            else
                % For more than 12 words, show every other word
                midPoints = (wordBoundaries(1:2:end-1) + wordBoundaries(2:2:end)) / 2;
                wordLabels = arrayfun(@(x) sprintf('Wd%d', x), 1:2:numWords, 'UniformOutput', false);
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

            % Determine overall maxval for proper scaling
            maxval = max(maxval_sentence, maxval_nonword);

            % Set enhanced y-axis limits
            ylim([-1.5, maxval + 0.5]);

            % Add significance markers
            if ~isempty(sigTimePoints)
                scatter(sigTimePoints, -1.25*ones(size(sigTimePoints)), 10, 'r', 'filled', 'LineWidth', 0.5);
            end

            % Add reference lines
            yline(0, 'k-.', 'LineWidth', 1, 'Alpha', 0.5);

            % Enhanced axis formatting
            xlabel('Word Position', 'FontSize', 8);
            ylabel('High Gamma Power', 'FontSize', 8);
            set(gca, 'FontSize', 10);
            grid on;
            grid minor;

            % Add trial count information
            text(0.02, 0.98, sprintf('Sentence trials: %d\nNonword trials: %d', ...
                size(trialData_sentence, 1), size(trialData_nonword, 1)), ...
                'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', 9, ...
                'BackgroundColor', 'white', 'EdgeColor', 'black');

            hold off;
        end

        % Enhanced legend
        legendAxes = axes('Position', [0.1, 0.01, 0.8, 0.08], 'Visible', 'off');
        hold(legendAxes, 'on');

        % Legend elements
        scatter(legendAxes, NaN, NaN, 10, 'r', 'filled');
        plot(legendAxes, NaN, NaN, 'Color', colors(1,:), 'LineWidth', 2);
        plot(legendAxes, NaN, NaN, 'Color', colors(2,:), 'LineWidth', 2);
        plot(legendAxes, NaN, NaN, 'k--', 'LineWidth', 1.5);
        plot(legendAxes, NaN, NaN, 'k-.', 'LineWidth', 1);

        legendLabels = {'Significant Sentence vs Nonword', ops.S_condition_flag, ops.N_condition_flag, ...
                       'Word Boundaries', 'Baseline'};
        legend(legendAxes, legendLabels, 'Orientation', 'horizontal', ...
               'Location', 'southoutside', 'FontSize', 8);

        hold(legendAxes, 'off');

        % Enhanced figure title
        sgtitle(sprintf('%s - Word Boundaries Analysis (UNIPOLAR) - Block %d/%d', ...
               obj.experiment, iChanBlock+1, totChanBlock), ...
               'FontSize', 8, 'FontWeight', 'bold');

        % Display the figure (or save it)
        if ops.do_plot
            set(f, 'Visible', 'on');
        end
    end

    % Process bipolar channels if available
    if isprop(obj, 's_vs_n_wordboundaries_concatenatedEpochsSentence_bip') && ...
       ~isempty(obj.s_vs_n_wordboundaries_concatenatedEpochsSentence_bip)

        concatenatedEpochsSentence_bip = obj.s_vs_n_wordboundaries_concatenatedEpochsSentence_bip;
        concatenatedEpochsNonword_bip = obj.s_vs_n_wordboundaries_concatenatedEpochsNonword_bip;
        pSig_bip = obj.s_vs_n_wordboundaries_pSig_bip;

        totChanBlock_bip = ceil(size(concatenatedEpochsSentence_bip, 1) / numChanBlock);

        % Similar plotting loop for bipolar channels (abbreviated for brevity)
        for iChanBlock = 0:(totChanBlock_bip-1)
            f = figure('Visible', 'off', 'Position', [100 100 1200 1400], 'Renderer', 'painters');

            for iChan = 1:min(numChanBlock, size(concatenatedEpochsSentence_bip, 1) - iChanBlock*numChanBlock)
                iChan2 = iChanBlock*numChanBlock + iChan;
                subplot(numChanBlock, 1, iChan);
                hold on;
    
                % Enhanced title with channel information
                title(sprintf('%s', obj.elec_ch_label{iChan2}), 'Interpreter', 'none', ...
                      'FontSize', 10, 'FontWeight', 'bold');
    
                % Process sentence data
                trialData_sentence = squeeze(concatenatedEpochsSentence_bip(iChan2, :, :));
                trialMean_sentence = nanmean(trialData_sentence, 1);
                trialSEM_sentence = nanstd(trialData_sentence, 0, 1) / sqrt(size(trialData_sentence, 1));
                maxval_sentence = max(trialMean_sentence + trialSEM_sentence);
    
                % Process nonword data
                trialData_nonword = squeeze(concatenatedEpochsNonword_bip(iChan2, :, :));
                trialMean_nonword = nanmean(trialData_nonword, 1);
                trialSEM_nonword = nanstd(trialData_nonword, 0, 1) / sqrt(size(trialData_nonword, 1));
                maxval_nonword = max(trialMean_nonword + trialSEM_nonword);
    
                % Plot sentence data with enhanced visualization
                plot(x, trialMean_sentence, 'Color', colors(1,:), 'LineWidth', 2);
                patch([x, fliplr(x)], [trialMean_sentence+trialSEM_sentence, ...
                      fliplr(trialMean_sentence-trialSEM_sentence)], ...
                      colors(1,:), 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    
                % Plot nonword data with enhanced visualization
                plot(x, trialMean_nonword, 'Color', colors(2,:), 'LineWidth', 2);
                patch([x, fliplr(x)], [trialMean_nonword+trialSEM_nonword, ...
                      fliplr(trialMean_nonword-trialSEM_nonword)], ...
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
                    wordLabels = arrayfun(@(x) sprintf('Wd%d', x), 1:numWords, 'UniformOutput', false);
                    set(gca, 'XTick', midPoints, 'XTickLabel', wordLabels);
                else
                    % For more than 12 words, show every other word
                    midPoints = (wordBoundaries(1:2:end-1) + wordBoundaries(2:2:end)) / 2;
                    wordLabels = arrayfun(@(x) sprintf('Wd%d', x), 1:2:numWords, 'UniformOutput', false);
                    set(gca, 'XTick', midPoints, 'XTickLabel', wordLabels);
                end
    
                % Plot significance markers with enhanced visualization
                sigTimePoints = [];
                for wordPos = 1:numWords
                    if wordPos <= size(pSig_bip, 2) && ~isempty(pSig_bip{iChan2, wordPos})
                        wordSigPoints = find(pSig_bip{iChan2, wordPos}.h_sig_05) + (wordPos-1)*timePointsPerWord;
                        sigTimePoints = [sigTimePoints, wordSigPoints];
                    end
                end
    
                % Determine overall maxval for proper scaling
                maxval = max(maxval_sentence, maxval_nonword);
    
                % Set enhanced y-axis limits
                ylim([-1.5, maxval + 0.5]);
    
                % Add significance markers
                if ~isempty(sigTimePoints)
                    scatter(sigTimePoints, -1.25*ones(size(sigTimePoints)), 10, 'r', 'filled', 'LineWidth', 0.5);
                end
    
                % Add reference lines
                yline(0, 'k-.', 'LineWidth', 1, 'Alpha', 0.5);
    
                % Enhanced axis formatting
                xlabel('Word Position', 'FontSize', 8);
                ylabel('High Gamma Power', 'FontSize', 8);
                set(gca, 'FontSize', 9);
                grid on;
                grid minor;
    
                % Add trial count information
                text(0.02, 0.98, sprintf('Sentence trials: %d\nNonword trials: %d', ...
                    size(trialData_sentence, 1), size(trialData_nonword, 1)), ...
                    'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', 9, ...
                    'BackgroundColor', 'white', 'EdgeColor', 'black');
    
                hold off;
            end
    
            % Enhanced legend
            legendAxes = axes('Position', [0.1, 0.01, 0.8, 0.08], 'Visible', 'off');
            hold(legendAxes, 'on');
    
            % Legend elements
            scatter(legendAxes, NaN, NaN, 30, 'r', 'filled');
            plot(legendAxes, NaN, NaN, 'Color', colors(1,:), 'LineWidth', 2);
            plot(legendAxes, NaN, NaN, 'Color', colors(2,:), 'LineWidth', 2);
            plot(legendAxes, NaN, NaN, 'k--', 'LineWidth', 1.5);
            plot(legendAxes, NaN, NaN, 'k-.', 'LineWidth', 1);
    
            legendLabels = {'Significant Sentence vs Nonword', ops.S_condition_flag, ops.N_condition_flag, ...
                           'Word Boundaries', 'Baseline'};
            legend(legendAxes, legendLabels, 'Orientation', 'horizontal', ...
                   'Location', 'southoutside', 'FontSize', 20);
    
            hold(legendAxes, 'off');

       

            sgtitle(sprintf('%s - Word Boundaries Analysis (BIPOLAR) - Block %d/%d', ...
                   obj.experiment, iChanBlock+1, totChanBlock_bip), ...
                   'FontSize', 20, 'FontWeight', 'bold');

            if ops.do_plot
                set(f, 'Visible', 'on');
            end
        end
    end

    % Display summary statistics
    fprintf(1, '> Unipolar analysis completed: %d channels, %d significant channels\n', ...
           size(concatenatedEpochsSentence, 1), length(obj.s_vs_n_wordboundaries_sigUnipolarChannels));

    if isfield(obj, 's_vs_n_wordboundaries_sigBipolarChannels')
        fprintf(1, '> Bipolar analysis completed: %d channels, %d significant channels\n', ...
               size(concatenatedEpochsSentence_bip, 1), length(obj.s_vs_n_wordboundaries_sigBipolarChannels));
    end
end
