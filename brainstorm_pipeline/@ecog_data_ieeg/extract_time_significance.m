function extract_time_significance(obj,args)
arguments
    obj ecog_data_ieeg
    args.probe_num = 1;
    args.baseTime = [-0.5 0];
    args.epochTime = [-0.5 6];
    args.numPerm = 1000;
    args.p_val = 0.05;
end
    baseTime = args.baseTime;
    epochTime = args.epochTime;

    [baseData,baseData_bip] = obj.extract_trial_epochs(epoch_tw = baseTime);
    [epochData,epochData_bip] = obj.extract_trial_epochs(epoch_tw = epochTime);

    baseDataExtend = extendTimeEpoch(baseData,size(epochData,3));

    

    parfor iChan = 1:size(baseDataExtend,1)
        aSig = squeeze(epochData(iChan,:,:));
        bSig = squeeze(baseDataExtend(iChan,:,:));
        pSigChan{iChan} = timePermCluster(aSig,bSig, pThresh=args.p_val);
    end
    obj.stats.time_series.pSigChan = pSigChan;
    if(~isempty(baseData_bip))
        baseData_bip_Extend = extendTimeEpoch(baseData_bip,size(epochData,3));
        parfor iChan = 1:size(baseData_bip_Extend,1)
            aSig = squeeze(epochData_bip(iChan,:,:));
            bSig = squeeze(baseData_bip_Extend(iChan,:,:));
            pSigChan_bip{iChan} = timePermCluster(aSig,bSig, pThresh=args.p_val);
        end
        obj.stats.time_series.pSigChan_bip = pSigChan_bip;
    end

    
    

end