%% 

%% writing ecog_data structure
subject = subjname;
experiment = taskname;
order = 'defaultSEEGorBOTHBroadBand';
save_filename = [ subject '_' experiment '_crunched_' order '.mat'];
save_path = [SaveDirectory filesep 'crunched' filesep];
currentDateTime = datetime('now');
formattedDateTime = datestr(currentDateTime, 'yyyymmdd_HHMM');
log_filename = ['../logs' filesep subject '_' experiment '_' order '_' formattedDateTime '.txt'];
d_files = nsxfileBank1;
ch_labels = ChanLabelBank1(GoodChannels1);
ch_nums = 1:length(ch_labels);
ch_type = {};
for iChan = 1:length(ch_labels)
    if(contains(ch_labels{iChan},ShankName1))
        ch_type{iChan} = 'seeg';
    else
        ch_type{iChan} = 'ground';
    end
end
ch_type = ch_type';


for_preproc = struct;
for_preproc.elec_data_raw = single(dataBank1(GoodChannels1,time2save));
for_preproc.event_table = event_table;
for_preproc.stitch_index_raw = 1;
for_preproc.sample_freq_raw = clinFS;
for_preproc.log_file_name = log_filename;
for_preproc.decimation_freq = clinFS/4;

obj = ecog_data(for_preproc,subject,experiment,save_filename,save_path,d_files,...
    DataDirectory,ch_labels,ch_nums,[],ch_type);
obj.preprocess_signal('order',order,'isPlotVisible',false,'doneVisualInspection',false);
obj.events_table = struct2table(obj.for_preproc.event_table);
obj.condition = extractfield(obj.for_preproc.event_table,'condition');
obj.session = extractfield(obj.for_preproc.event_table,'session');

obj.trial_timing = make_trials_swm(obj.events_table,clinFS);


if(~isfolder(save_path))
    mkdir(save_path)
end

save([save_path filesep save_filename],'obj','-v7.3');

%% Additional analysis by Kumar


% Extracting High Gamma Signals
obj.extract_high_gamma('doNapLabFilterExtraction',true);
% Downsampling to 200Hz
obj.downsample_signal('decimationFreq',200);
% Finding significant channels using windowed averaging permutation test
obj.extract_significant_channel();
% Finding significant channels using cluster correction based time-series
% permutation test - optional & takes a lot of time
obj.extract_time_significance();
% Extracting baseline normalization metrics
obj.extract_normalization_metrics();
% Z-scoring using baseline normalization metrics
obj.normalize_signal("normtype",'z-score');
% Extracting trial epoch
[epochData,epochData_bip] = obj.extract_trial_epochs(epoch_tw = [-0.5 6]);
% 
figure; plot(squeeze(mean(epochData_bip,2))')
sig_ch = obj.stats.sig_hg_channel.h_bip_fdr_05;
figure; plot(squeeze(mean(epochData_bip(sig_ch,:,:),2))')

sig_ch_id = find(sig_ch);

hTrials = get_cond_id(obj,'H');
eTrials = get_cond_id(obj,'E');
for iChan = 1:length(sig_ch_id)
    figure;
    plot(squeeze(nanmean(epochData_bip(sig_ch_id(iChan),eTrials,:),2))');
    hold on;
    plot(squeeze(nanmean(epochData_bip(sig_ch_id(iChan),hTrials,:),2))');
end
