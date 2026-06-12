function [Xh_out,X_fft_h_out]=hilbert_transform(obj,X,rate,varargin)
    % https://github.com/BouchardLab/ecogVIS/blob/master/ecogvis/signal_processing/hilbert_transform.py
    % Based on : Computing the Discrete-Time Analytic Signal via FFTS.
    % Lawrence  Marple,  Jr. IEEE 1999
    %
    % Arguments
    % ---------
    % X : ndarray (n_channels, n_time)
    %   Input data, dimensions
    % rate : float
    %   Number of samples per second.
    % filters : filter or cell of filters (optional)
    %   One or more bandpass filters
    % 
    % Returns
    % -------
    % Xh : ndarray, complex
    %   Bandpassed analytic signal
    %

    switch nargin 
    case 4
        filters=varargin{1};
        phase='None';
        X_fft_h='None';
    case 5
        filters=varargin{1};
        phase=varargin{2};
        X_fft_h='None';
    case 6
        filters=varargin{1};
        phase=varargin{2};
        X_fft_h=varargin{3};
    otherwise 
        filters='None';
        phase='None';
        X_fft_h='None';        
    end 
        
    if ~iscell(filters)
        filters = {filters};
    end 

    time = size(X,2);
    d=1./rate;
    freq=[0:ceil(time/2-1), ceil(-(time)/2):1:-1]/(d*time);
    Xh=cell(length(filters),1);
    %Xh = zeros((len(filters),) + X.shape, dtype=np.complex)

    % compute the one sided analytic signal in frequency domain 
    if strcmp(X_fft_h, 'None')
        % Heavyside filter
        h = zeros(1,length(freq));
        h(freq > 0) = 2.;
        h(1) = 1.;
        % X_fft_h is the single sided spectrum of X
        X_fft_h = fft(X) .* h;

        if ~strcmp(phase,'None')
            X_fft_h=X_fft.*phase;
        end
    end
        
    for ii=1:size(filters,2)
        f=filters{ii};
        if strcmp(f,'None')
            Xh{ii} = ifft(X_fft_h);
        else
            % filter the signal in different frequency bands 
            f = f / max(f);
                
            Xh{ii} = ifft(X_fft_h .* f);
        end 
    end

    if size(Xh,1) == 1
        Xh_out= Xh{1};
        X_fft_h_out=X_fft_h;
    else 
        Xh_out= Xh;
        X_fft_h_out=X_fft_h;
    end 
    
end