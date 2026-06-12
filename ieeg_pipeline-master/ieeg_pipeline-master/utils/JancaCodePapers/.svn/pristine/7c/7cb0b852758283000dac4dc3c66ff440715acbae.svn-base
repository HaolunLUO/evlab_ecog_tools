%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% JB Eichenlaub - Last modification 11.03.2016
% use spike_detector_hilbert_v21 script from Janca et al. 

function [spikesDetected] = automaticSpikeDetection_UsingJancaMethod(d, fs, settings)

% v21. Janca et al., 2015, Brain Topogr - "Detection of Interictal Epileptiform ..."
[out, discharges, d_decim, envelope, background, envelope_pdf] = spike_detector_hilbert_v21(d, fs, settings);

% Extract Data
spikesDetected.out          = out;
spikesDetected.discharges   = discharges;
spikesDetected.d_decim      = d_decim;
spikesDetected.envelope     = envelope;
spikesDetected.background   = background;
spikesDetected.envelope_pdf = envelope_pdf;

clearvars -except spikesDetected
