function ieegPad = extendTimeEpoch(ieegData, sigLen)
    warning off;
    ieegPad = [];
    % Validate input data dimensions
    if ndims(ieegData) < 3
        error('ieegData must have at least three dimensions.');
    end

    % Generate padTimeArray safely
    padTimeArray = 50:50:size(ieegData,3);

    if isempty(padTimeArray)
        padTimeArray = 5:5:size(ieegData,3);
    end

    if isempty(padTimeArray)
        error('padTimeArray is empty.');
    end

    %try
        for iChan = 1:size(ieegData,1)
            ieegChan = squeeze(ieegData(iChan, :, :));
            ieegChanPad = [];

            for iTrial = 1:size(ieegChan, 1)
                validPadChoose = false;
                availablePadTimes = padTimeArray; % copy to modify
                sampleSize = 0;

                while ~validPadChoose && ~isempty(availablePadTimes)
                    padChoose = randsample(availablePadTimes,1);
                    time2pad = sigLen / padChoose;

                    % if mod(sigLen, padChoose) ~= 0
                    %     warning('sigLen is not divisible by padChoose. Adjusting time2pad.');
                    %     time2pad = floor(time2pad);
                    % end

                    selectTrials = setdiff(1:size(ieegChan, 1), iTrial);
                    if isempty(selectTrials)
                        error('No trials available for padding.');
                    end

                    sampleSize = ceil(time2pad) - 1;

                    if sampleSize <= numel(selectTrials)
                        validPadChoose = true;
                    else
                        % Remove this padChoose and try another
                        availablePadTimes(availablePadTimes == padChoose) = [];
                    end
                end

                % If no valid padChoose found, use best possible
                if ~validPadChoose
                    sampleSize = min(sampleSize, numel(selectTrials));
                end

                if sampleSize > 0
                    randTrials = datasample(selectTrials, sampleSize, 'Replace', false);
                else
                    randTrials = [];
                end

                if ~isempty(randTrials)
                    trials2join = ieegChan(randTrials, 1:padChoose)';

                    if ceil(time2pad) == time2pad
                        ieegChanPad(iTrial, :) = [ieegChan(iTrial, 1:padChoose) trials2join(:)'];
                    else
                        joinTrials = trials2join(:)';
                        ieegChanPad(iTrial, :) = [ieegChan(iTrial, 1:padChoose) joinTrials(1:sigLen - padChoose)];
                    end
                else
                    % Pad with random noise of mean 0 and std 0.25
                    padLength = sigLen - padChoose;
                    randomPad = randn(1, padLength) * 0.25;
                    ieegChanPad(iTrial, :) = [ieegChan(iTrial, 1:padChoose) randomPad];
                end
            end

            ieegPad(iChan,:,:) = ieegChanPad;
        end

    % catch ME
    %     fprintf('Error occurred: %s\n', ME.message);
    % end
end
