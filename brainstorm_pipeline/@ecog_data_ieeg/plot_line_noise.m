function f=plot_line_noise(obj,noise_before,noise_after)
    % 

    fprintf(1, '\n> Plotting 60Hz noise power ...\n');

    close all
    f = figure;
    set(gcf,'position',[30,30,2300,900]);
    c= [0.4660 0.6740 0.1880];

    x = find(~obj.elec_ch_valid);
    idxs = ~obj.elec_ch_valid;

    if(isempty(x))
        fprintf(1, '\n> No Invalid channels detected')
    else

    % --- ALL NOISE BEFORE ---
    currsub = subplot(2,2,1);       
    stem(noise_before,'filled'); 
    axis tight; hold on;
    if(isscalar(x))
        [nPoints, nTime] = size(noise_before(idxs, :));
        x_expanded = repmat(x, 1, nTime);
    else
        x_expanded = x;
    end
    stem(x_expanded,noise_before(idxs,:),'filled','Color','k')    
    legend({'55Hz noise','60Hz noise','65Hz noise','MARKED NOISY'},'Location','best','FontSize',16,'Box','off');
    ylabel('Noise (uV)','FontSize',18);
    title('BEFORE NOTCH FILTERING','FontSize',22);
    obj.update_position(currsub);

    % --- NOISE RATIO BEFORE ---
    currsub = subplot(2,2,3);
    stem(noise_before(:,2)./mean(noise_before(:,[1,3]),2),'filled','Color',c); 
    axis tight; hold on;
    stem(x,noise_before(idxs,2)./mean(noise_before(idxs,[1,3]),2),'filled','Color','k')
    xlabel('Channel #','FontSize',18);
    ylabel('60Hz noise / mean 55Hz+65Hz noise','FontSize',18)
    obj.update_position(currsub);

    % --- ALL NOISE AFTER ---
    currsub = subplot(2,2,2);
    stem(noise_after,'filled'); 
    axis tight; hold on;
    stem(x_expanded,noise_after(idxs,:),'filled','Color','k')
    title('AFTER NOTCH FILTERING','FontSize',22)
    obj.update_position(currsub);

    % --- NOISE RATIO AFTER ---
    currsub = subplot(2,2,4);
    stem(noise_after(:,2)./mean(noise_after(:,[1,3]),2),'filled','Color',c); 
    axis tight; hold on;
    stem(x,noise_after(idxs,2)./mean(noise_after(idxs,[1,3]),2),'filled','Color','k')
    xlabel('Channel #','FontSize',18)
    obj.update_position(currsub);

    % --- SAVE ---
    PATH = [obj.crunched_file_path 'plots/line_noise/'];
    if(~isdir(PATH))
        mkdir(PATH);
    end
    filename = split(obj.for_preproc.log_file_name,'/');
    filename = split(filename{end},'.');
    filename = [filename{1} '_line_noise_invalid_channels.png'];
    saveas(gcf,strcat(PATH,filename));
    set(0, 'CurrentFigure', f);
    end


    x = find(obj.elec_ch_valid);
    idxs = obj.elec_ch_valid;

    if(isempty(x))
        fprintf(1, '\n> No valid channels detected')
    else

    % --- ALL NOISE BEFORE ---
    currsub = subplot(2,2,1);       
    stem(noise_before,'filled'); 
    axis tight; hold on;
    stem(x,noise_before(idxs,:),'filled','Color','k')
    legend({'55Hz noise','60Hz noise','65Hz noise','MARKED NOISY'},'Location','best','FontSize',16,'Box','off');
    ylabel('Noise (uV)','FontSize',18);
    title('BEFORE NOTCH FILTERING','FontSize',22);
    obj.update_position(currsub);

    % --- NOISE RATIO BEFORE ---
    currsub = subplot(2,2,3);
    stem(noise_before(:,2)./mean(noise_before(:,[1,3]),2),'filled','Color',c); 
    axis tight; hold on;
    stem(x,noise_before(idxs,2)./mean(noise_before(idxs,[1,3]),2),'filled','Color','k')
    xlabel('Channel #','FontSize',18);
    ylabel('60Hz noise / mean 55Hz+65Hz noise','FontSize',18)
    obj.update_position(currsub);

    % --- ALL NOISE AFTER ---
    currsub = subplot(2,2,2);
    stem(noise_after,'filled'); 
    axis tight; hold on;
    stem(x,noise_after(idxs,:),'filled','Color','k')
    title('AFTER NOTCH FILTERING','FontSize',22)
    obj.update_position(currsub);

    % --- NOISE RATIO AFTER ---
    currsub = subplot(2,2,4);
    stem(noise_after(:,2)./mean(noise_after(:,[1,3]),2),'filled','Color',c); 
    axis tight; hold on;
    stem(x,noise_after(idxs,2)./mean(noise_after(idxs,[1,3]),2),'filled','Color','k')
    xlabel('Channel #','FontSize',18)
    obj.update_position(currsub);

    % --- SAVE ---
    PATH = [obj.crunched_file_path 'plots/line_noise/'];
    if(~isdir(PATH))
        mkdir(PATH);
    end
    filename = split(obj.for_preproc.log_file_name,'/');
    filename = split(filename{end},'.');
    filename = [filename{1} '_line_noise_valid_channels.png'];
    saveas(gcf,strcat(PATH,filename));
    set(0, 'CurrentFigure', f);
    end

end