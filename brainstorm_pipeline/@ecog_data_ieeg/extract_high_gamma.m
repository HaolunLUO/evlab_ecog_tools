function extract_high_gamma(obj,ops)
    % Extracts high gamma signal using one of two methods: 
    %   (1) Chang Lab method (gaussian filters)
    %           hilbert_transform.m is a special function 
    %           See 'HELPER SIGNAL PREPROCESSING FUNCTIONS'
    %   (2) Standard bandpass filter
    %           hilbert.m is a standard MATLAB function
    % 
    % TODO - descriptions of above methods

    arguments
        obj ecog_data_ieeg
        ops.doGaussianFilterExtraction = false; % similar to Chang lab filters but uses slightly varying parameters
        ops.doBandpassExtraction = false;
        ops.doNapLabFilterExtraction = false; % similar to Chang lab filters but uses modules from NAPLAB library
    end

    % p = inputParser();
    % addParameter(p,'doGaussianFilterExtraction',false);
    % addParameter(p,'doBandpassExtraction',false);
    % addParameter(p,'doNapLabFilterExtraction', false); % similar to Chang lab filters but uses slightly varying parameters
    % parse(p, varargin{:});
    % ops = p.Results;

    % ensure that at least one, and only one, extraction method is specified
    if ops.doGaussianFilterExtraction == ops.doBandpassExtraction
        if ops.doBandpassExtraction == ops.doNapLabFilterExtraction
            error('Must specify one - and only one - method for extracting the high gamme signal')
    
        end
    end

    signal = double(obj.elec_data');
    if ~isempty(obj.bip_elec_data)
        signal_bipolar = obj.bip_elec_data';
    end

    % define filter params if it hasn't already been done
    if ~isfield(obj.for_preproc,'gaussian') || ~isfield(obj.for_preproc,'bandpass')
        obj.define_parameters()
    end


    for k=1:length(obj.stitch_index) % number of separate data files with signal
        fprintf(1, '\n> Extracting high gamma signal from file %d of %d ... \n',k,length(obj.stitch_index));

        if k == length(obj.stitch_index) % signal for file stops at end of matrix
            stop = size(signal,1);
        else % signal for file stops before stitch index of next file
            stop = obj.stitch_index(k+1)-1;
        end

        signal_ = signal(obj.stitch_index(k):stop,:);


        % ------------------------------------
        % EXTRACTION USING GUASSIAN FILTERING
        % ------------------------------------
        if ops.doGaussianFilterExtraction

            cfs = obj.for_preproc.gaussian.cfs;
            sds = obj.for_preproc.gaussian.sds;

            % define gaussian filters
            % filter_bank = cell(numel(cfs),1);
            filter_bank={};
            for s=1:length(cfs)
                filter_bank{s} = obj.gaussian_filter(transpose(signal_(:,1)),obj.sample_freq,cfs(s),sds(s));
            end
            obj.for_preproc.gaussian.filter_banks = filter_bank;
            % size(filter_bank)

            % --- UNIPOLAR ---
            fprintf(1, '\n>> Extracting unipolar high gamma envelope based on gaussian filtering \n');
            fprintf(1,'[');

            signal_hilbert = nan*signal_;

            for kk=1:size(signal_,2)
                signal_hilbert_all = cell2mat(cellfun(@abs,obj.hilbert_transform(double(transpose(signal_(:,kk))),obj.sample_freq,filter_bank),'UniformOutput',false));
                % size(signal_hilbert_all)
                signal_hilbert(:,kk) = transpose(mean(signal_hilbert_all,1));

                fprintf(1,'.');

            end 

            fprintf(1,'] done\n');

            signal(obj.stitch_index(k):stop,:) = signal_hilbert;


            % --- BIPOLAR ---
            if ~isempty(obj.bip_elec_data)
                fprintf(1, '\n>> Extracting bipolar high gamma envelope based on gaussian filtering \n');
                fprintf(1,'[');

                signal_bipolar_ = signal_bipolar(obj.stitch_index(k):stop,:);
                signal_hilbert_bipolar = nan*signal_bipolar_;

                for kk=1:size(signal_bipolar_,2)
                    signal_hilbert_bipolar_all = cell2mat(cellfun(@abs,obj.hilbert_transform(double(transpose(signal_bipolar_(:,kk))),obj.sample_freq,filter_bank),'UniformOutput',false));
                    
                    signal_hilbert_bipolar(:,kk) = transpose(mean(signal_hilbert_bipolar_all,1));

                    fprintf(1,'.');

                end

                fprintf(1,'] done\n')

                signal_bipolar(obj.stitch_index(k):stop,:) = signal_hilbert_bipolar;
                
            end

        % ---------------------------------------------
        % EXTRACTION USING NAPLAB GAUSSIAN FILTERING
        % ---------------------------------------------
        elseif ops.doNapLabFilterExtraction
            freqRange = [obj.for_preproc.filter_params.bandpass.f_gamma_low obj.for_preproc.filter_params.bandpass.f_gamma_high];
             % --- UNIPOLAR ---
            fprintf(1, '\n>> Extracting unipolar high gamma envelope based on NAPLAB gaussian filtering \n');
            fprintf(1,'[');

            

            [dh,cfs,sigma_fs] = CUprocessingHilbertTransform_filterbankGUI(signal_', obj.sample_freq, freqRange);
            signal_hilbert = mean(abs(dh),3); 
            
            fprintf(1,'] done\n');

            signal(obj.stitch_index(k):stop,:) = signal_hilbert';


            % --- BIPOLAR ---
            if ~isempty(obj.bip_elec_data)
                fprintf(1, '\n>> Extracting bipolar high gamma envelope based on NAPLAB gaussian filtering \n');
                fprintf(1,'[');

                signal_bipolar_ = signal_bipolar(obj.stitch_index(k):stop,:);
                

                [dh,cfs,sigma_fs] = CUprocessingHilbertTransform_filterbankGUI(signal_bipolar_', obj.sample_freq, freqRange);
                signal_hilbert_bipolar = mean(abs(dh),3);

                

                signal_bipolar(obj.stitch_index(k):stop,:) = signal_hilbert_bipolar';
                
            end

        % ---------------------------------------------
        % EXTRACTION USING STANDARD BANDPASS FILTERING
        % ---------------------------------------------
        elseif ops.doBandpassExtraction

            B = obj.for_preproc.bandpass.B;
            A = obj.for_preproc.bandpass.A;

            % --- UNIPOLAR ---
            fprintf(1, '\n>> Extracting unipolar high gamma envelope based on bandpass filtering \n');

            % apply filter
            signal_hilbert = filtfilt(B,A,double(signal_));

            % measure envelope
            signal_hilbert = abs(hilbert(signal_hilbert));

            % truncate
            signal_hilbert(signal_hilbert < 0) = 0;


            % --- BIPOLAR ---
            if ~isempty(obj.bip_elec_data)
                fprintf(1, '\n>> Extracting bipolar high gamma envelope based on bandpass filtering \n');

                signal_hilbert_bipolar = filtfilt(B,A,double(signal_bipolar_));
                signal_hilbert = abs(hilbert(signal_hilbert));
                signal_hilbert(signal_hilbert < 0) = 0;

                signal_bipolar(obj.stitch_index(k):stop,:) = signal_hilbert_bipolar;
            end
        
        end

    end

    obj.elec_data = signal';

    if ~isempty(obj.bip_elec_data)
        obj.bip_elec_data = signal_bipolar';
    end

end



function [filteredData,cfs,sigma_fs,hilbdata]=CUprocessingHilbertTransform_filterbankGUI(d,Fs,freqRange)

% This function is used in EcogExtractHighGamma.m
%{
PURPOSE: Perform Hilbert transform

INPUTS: ecog data structure
        Sampling frequency
        Optional: frequency range for window (2 element array with low frequency first)-- if no input, go to default range

OUTPUT: filtered data structure
%}

%********CHANGE, allow for multiband freqRange*************
%{
if nargin==3
    freqH=freqRange(2);
    freqL=freqRange(1);
else
    freqH=150;
    freqL=70;
end

max_freq=Fs/2;
%}
% Neural Acoustic Processing Lab, 
% Columbia University, naplab.ee.columbia.edu


%%%%%%%%%%%%%%%CREATE FILTER BANK
a=[log10(.39); .5];

frange=freqRange;
f0=0.018;
octspace=1/7;%usually 1/7
minf=frange(1);
maxf=frange(2);
maxfo=log2(maxf/f0);
cfs=f0;
sigma_f=10^(a(1)+a(2)*log10(cfs(end)));

while log2(cfs(end)/f0)<maxfo
    cfo=log2(cfs(end)/f0);
    cfo=cfo+octspace;
    if cfs(end)<4,
        cfs=[cfs cfs(end)+sigma_f]; %switches to log spacing at 4 Hz
    else cfs=[cfs f0*(2^(cfo))];
    end
    sigma_f=10^(a(1)+a(2)*log10(cfs(end)));
end

cfs=cfs(find(cfs>=minf & cfs<=maxf));
npbs=length(cfs);
sigma_fs=(10.^([ones(length(cfs),1) log10(cfs')]*a))';
badfs=[find(cfs>340 & cfs<480) find(cfs>720 & cfs<890)];
sigma_fs=sigma_fs(setdiff(1:npbs,badfs));
cfs_all=cfs;
cfs=cfs(setdiff(1:npbs,badfs));
npbs=length(cfs);
sds=sigma_fs.*sqrt(2);

T=size(d,2);
freqs=(0:floor(T/2)).*(Fs/T); nfreqs=length(freqs);
h = zeros(1,T);
if 2*fix(T/2)==T %if T is even
    h([1 T/2+1]) = 1;
    h(2:T/2) = 2;
else
    h(1) = 1; h(2:(T+1)/2) = 2;
end

%CHANGE: vectorize across channels, take out loop*******************
%x=fft(ecog.data,nfft,2);
filteredData = zeros(size(d,1),T,npbs);
fprintf(1,'[');
for c=1:size(d,1)
    adat=fft(d(c,:),T);
    for f=1:npbs
        H = zeros(1,T);
        k = freqs-cfs(f);
        H(1:nfreqs) = exp((-0.5).*((k./sds(f)).^2));
        H(nfreqs+1:end)=fliplr(H(2:ceil(T/2)));
        H(1)=0;
        hilbdata=ifft(adat(end,:).*(H.*h),T);
        envData=abs(hilbdata);
        %phaseData=angle(hilbdata);
        filteredData(c,:,f)=hilbdata;
        %phaseInfo.data(c,:,f)=phaseData;
    end
    fprintf(1,'.');
end
fprintf(1,'] done\n')
filteredData = abs(filteredData);

end