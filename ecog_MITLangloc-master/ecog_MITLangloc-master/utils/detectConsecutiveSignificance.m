function [has_consecutive_sig, num_consecutive, first_consecutive_word] = ...
    detectConsecutiveSignificance(pSig_channel, nWords, min_consecutive)
    % Detects if channel has significance in min_consecutive (default 3) consecutive words
    %
    % Inputs:
    %   pSig_channel     - Cell array of significance results for one channel across words
    %   nWords           - Total number of words
    %   min_consecutive  - Minimum number of consecutive significant words (default: 3)
    %
    % Outputs:
    %   has_consecutive_sig   - Logical: true if consecutive significant words found
    %   num_consecutive       - Number of consecutive significant words (max sequence)
    %   first_consecutive_word - First word index of the longest consecutive sequence
    
    if nargin < 3
        min_consecutive = 3;  % Default: require 3 consecutive words
    end
    
    % Initialize flags for each word position
    sig_per_word = false(nWords, 1);
    
    % Check each word position for significance
    for wordPos = 1:nWords
        if ~isempty(pSig_channel{wordPos})
            % Check if any timepoint in this word shows significance
            if any(pSig_channel{wordPos}.h_sig_05 == 1)
                sig_per_word(wordPos) = true;
            end
        end
    end
    
    % Find consecutive sequences of significant words
    has_consecutive_sig = false;
    num_consecutive = 0;
    first_consecutive_word = 0;
    max_consecutive = 0;
    max_consecutive_start = 0;
    
    % Scan for consecutive significant words
    current_consecutive = 0;
    current_start = 0;
    
    for wordPos = 1:nWords
        if sig_per_word(wordPos)
            if current_consecutive == 0
                % Start of a new consecutive sequence
                current_start = wordPos;
            end
            current_consecutive = current_consecutive + 1;
        else
            % End of consecutive sequence
            if current_consecutive > max_consecutive
                max_consecutive = current_consecutive;
                max_consecutive_start = current_start;
            end
            current_consecutive = 0;
        end
    end
    
    % Check the last sequence
    if current_consecutive > max_consecutive
        max_consecutive = current_consecutive;
        max_consecutive_start = current_start;
    end
    
    % Determine if threshold is met
    if max_consecutive >= min_consecutive
        has_consecutive_sig = true;
        num_consecutive = max_consecutive;
        first_consecutive_word = max_consecutive_start;
    end
end