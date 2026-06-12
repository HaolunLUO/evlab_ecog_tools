function plot_channels(obj,varargin)
    % TODO - description

    

    p = inputParser();
    addRequired(p,'signal');
    addRequired(p,'channel_labels');
    addRequired(p,'clean_channels');
    addRequired(p,'valid_channels');
    addParameter(p,'stitch_index',[1]);
    addParameter(p,'t_len',60); % size of viewing window in seconds
    addParameter(p,'sample_freq',1200);
    addParameter(p,'downsample',false);
    addParameter(p,'decimation_freq',300)
    addParameter(p,'plotIEDs',false);
    addParameter(p,'chanIEDs',[]);
    addParameter(p,'posIEDs',[]);
    addParameter(p,'save',false); 
    parse(p, varargin{:});
    ops = p.Results;

    curr_sample_freq = ops.sample_freq;
    D = ops.signal;


    % ------------------------------
    % FORMAT & NORMALIZE SIGNAL
    % ------------------------------
    x_norm_cell = [];
    for k=1:length(ops.stitch_index) % number of separate data files with signal

        if k == length(ops.stitch_index) % signal for file stops at end of matrix
            stop = size(D,1);
        else % signal for file stops before stitch index of next file
            stop = ops.stitch_index(k+1)-1;
        end 

        D_ = D(ops.stitch_index(k):stop,:);

        % formatting
        x_cell = mat2cell(D_',ones(1,size(D_,2))); % make signal from each electrode a cell

        % downsampling
        if ops.downsample
            decimation_factor = ops.sample_freq/ops.decimation_freq;
            x_cell = cellfun(@(x) downsample(x,decimation_factor),x_cell,'uni',false);
            curr_sample_freq = ops.decimation_freq;
        end

        % normalizing
        min_max = mean(cell2mat(cellfun(@(y) prctile(y,[3 97]),x_cell,'UniformOutput',false))); % mean 5th and 95th percentiles of ALL electrodes
        x_norm_cell_ = cellfun(@(x) (x-min_max(1))./(min_max(2)-min_max(1)),x_cell,'UniformOutput',false); % normalize
        x_norm_cell_ = arrayfun(@(x) x_norm_cell_{x}*ops.valid_channels(x),1:size(x_norm_cell_,1),'uni',false)'; % set signal of noisy channels to 0
        
        if length(ops.stitch_index) > 1
            x_norm_cell = [x_norm_cell, x_norm_cell_];
        else 
            x_norm_cell = x_norm_cell_;
        end

    end

    % combine separate x_norm_cell columns
    if length(ops.stitch_index) > 1
        for k=1:length(obj.stitch_index)-1
            x_norm_cell(:,k+1) = arrayfun(@(x) {[x_norm_cell{x,k}, x_norm_cell{x,k+1}]},[1:size(x_norm_cell,1)])';
        end
        x_norm_cell = x_norm_cell(:,length(ops.stitch_index));
    end
    assert(size(x_norm_cell,2)==1,'x_norm_cell not in the correct format');


    % ------------
    % PLOT SIGNAL
    % ------------
    t_length = ops.t_len*curr_sample_freq;

    col_inf=inferno(floor(.8*size(x_norm_cell,1)));
    col_vir=viridis(floor(.8*size(x_norm_cell,1)));
    colors=[col_vir(1:floor(size(x_norm_cell,1)/2),:);col_inf(1:(floor(size(x_norm_cell,1)/2)+1),:)];

    close all
    figure(1);
    clf;
    set(gcf,'position',[31,1,1700,900]); % 1713, 1010
    ax = axes('position',[.05,.1,.93,.88]);
    hold on
    time_stamps = [1:size(x_norm_cell{1},2)]/curr_sample_freq;
    hold on
    H = arrayfun(@(x) plot(time_stamps,x_norm_cell{x}+x,'color',colors(x,:),'tag',sprintf('ch %d, tag %s',x,ops.channel_labels{x})),[1:size(x_norm_cell,1)]);
       

    if ops.plotIEDs
        spike_chan = ops.chanIEDs;
        [spike_chan_sort,sort_idx] = sort(spike_chan);
        spike_times_sort = ops.posIEDs(sort_idx);
        
        yval = cell2mat(arrayfun(@(x) x_norm_cell{spike_chan_sort(x)}(floor(spike_times_sort(x)*200))+spike_chan_sort(x),1:length(spike_chan_sort),'uni',false));
        H1 = scatter(spike_times_sort,yval,50,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 .7 .7],'LineWidth',2,'MarkerFaceAlpha',.5);
    end
        
    set(gcf,'doublebuffer','on');
    set(ax,'ytick',[1:size(x_norm_cell,1)]);
    set(ax,'yticklabel','');
    arrayfun(@(x) text(0,x,num2str(x),'Color',colors(x,:),'HorizontalAlignment','right','VerticalAlignment','middle'),ops.clean_channels);
       
    set(ax,'ylim',[0,size(x_norm_cell,1)+4]);
    ax.XAxis.TickLength = [0.005,0.01];
    ax.YAxis.TickLength = [0.005,0.01];
    set(ax,'xlim',[0 ops.t_len]);
    pos = get(ax,'position');
    Newpos = [pos(1) pos(2)-0.1 pos(3) 0.05];
    xmax=max(time_stamps);
    S = ['set(gca,''xlim'',get(gcbo,''value'')+[0 ' num2str(ops.t_len) '])'];
        
    h = uicontrol('style','slider','units','normalized','position',Newpos,'callback',S,'min',0,'max',xmax-ops.t_len);
    datacursormode on
        
    waitfor(findobj('type','figure','number',1));

end