function SBA_SPC() 

    clc
    close all
        
    %% Main data directory
    global Dir
    Dir = '/Users/sudipvhaduri/Desktop/IEEE TBD/share/data/';
        
    %% list of subject IDs
    global subjList
    subjList=1:5;
    
    %% Set clustering parameters 
    global Tmax Tmin Dmax
    Tmax = 10; % minutes
    Tmin = 10; % minutes
    Dmax = 250; % meters
        
    %% Pairing sp-clusters with battery and sleep info 
    algo_type = {'SSPC','BSPC','SBSPC'}; 
%     parfor i=1:size(algo_type,2)
    for i=1:size(algo_type,2)
        disp([algo_type{1,i} ' algorithm is running ...'])
        combine_sps_using_sleep_battery_info(algo_type{1,i});
    end
end
function combine_sps_using_sleep_battery_info(algo_type)
    global Dir
    global subjList
    global Tmax Tmin Dmax
    
    outDir = [Dir 'outDir/SubjectLevel_sp_median_pairs_Tmin' num2str(Tmin) 'minutes_Dmax' num2str(Dmax) 'meters_BatteryExtensionsAt99BL/'];
    if(~exist(outDir,'dir')); mkdir(outDir); end
    Dir_sp = [Dir 'inDir/SubjectLevel_sp_median_Tmin' num2str(Tmin) 'minutes_Dmax' num2str(Dmax) 'meters_modifiedLocations/'];
    Dir_BS = [Dir 'inDir/SubjectLevel_plugUnplug_Seg_with_ovChg_Mkr_SlpMeta_overlapping_Extended_with_Place_Labels/'];
    Dir_SleepMeta = [Dir 'inDir/SubjectLevel_SleepMeta_Sleep_Segments_with_Place_Labels/'];
       
    outDir_SeqMerge = [Dir 'outDir/SubjectLevel_sp_median_pairs_Tmin' num2str(Tmin) 'minutes_Dmax' num2str(Dmax) 'meters_BatteryExtensionsAt99BL_SequenceMerging_' algo_type '/'];
    if(~exist(outDir_SeqMerge,'dir')) mkdir(outDir_SeqMerge); end;
    
    noise_affected_clustering_all_subj = [];
    BS_available_all_subj = [];
    Slp_available_all_subj = [];
    parfor i=1:length(subjList)
%     for i=1:length(subjList)
        pid = subjList(i) 
        
        %% Read all the sp-clusters of an individual 
        FILE_Name = ['sp_pid_' num2str(pid)];
        FILEs = [Dir_sp FILE_Name '.csv'];
        if(~exist(FILEs, 'file')); continue; end
        fid2 = fopen(FILEs);
        % sp#,c_lat_med,c_long_med,st_ts,end_ts
        B = textscan(fid2,'%d%f%f%f%f','delimiter',',','headerlines',1);
        fclose(fid2);
        sp_no=B{1,1};
        if(length(sp_no)<2); continue; end
        [c_lat_med,c_long_med]=B{1,2:3}; 
        lat_long_sp = [c_lat_med,c_long_med];
        [st_ts_sp,end_ts_sp]=B{1,4:5}; ts_sp = [st_ts_sp end_ts_sp]; % ts are in sec
        ts_sp = ts_sp.*1000; % in mili-sec
        
        %% Read all SleepMeta Sleep Sessions with place labels of an individual
        File = [Dir_SleepMeta 'NetHealth_' num2str(pid) '_Sleep.csv'];
        ts_s = [];
        if(exist(File,'file')) 
            fid = fopen(File);
            % ts_start(mili-sec),ts_end(mili-sec)
            B = textscan(fid,'%f%f','delimiter',',','headerlines',1);
            fclose(fid);
            [ts_st_s,ts_end_s]=B{1,1:2}; % mili-sec 
            ts_s = [ts_st_s ts_end_s];
        end
        
        %% Read all Battery Charging Sessions with place labels of an individual
        File = [Dir_BS 'NetHealth_' num2str(pid) '_plugUnplug_ext' '.csv'];
        ts_b = []; 
        if(exist(File,'file'))
            fid = fopen(File);
            % st_ts,end_ts,extended_end_ts,BL(percentage during extension)
            B = textscan(fid,'%f%f%f%d','delimiter',',','headerlines',1);
            fclose(fid);
            [ts_st_b,ts_end_b,ts_end_ext_b]=B{1,1:3}; 
            BL=B{1,4};
            ts_b = [ts_st_b ts_end_b];
            ind_temp = find(ts_end_ext_b>0 & BL==99); % consider those extensions that observe callback @99% battery level
            ts_b(ind_temp,2) = ts_end_ext_b(ind_temp);
        end
        
        %% sp cluster pairs with battery and sleep info
        OutFile = [outDir 'sp_pairs_pid_' num2str(pid) '.csv'];
        fileID = fopen(OutFile,'w+');
        fprintf(fileID,'sp#1,sp#2,ds(meter-bet.medianCentroids),dt(minutes-bet.spBoundaries),on_off(both sp on-campus?),case#(using battery),uncovered_gap_bet_sp(in minutes|using battery),case#(using sleep), uncovered_gap_bet_sp(in minutes|using sleep),timeExtension(using battery),timeExtension(using sleep), merged(-1=NO INFO|0=NOT Merged|1=BOTH battery and sleep merged|2=battery merged|3=sleep merged|4=battery and sleep together merged),timeExtension(2values=merged|4values=extendedButNotMerged|empty=NotExtended),actualSpPairTime\n');
        
        %% find inter-sp gap (both time and distance)
        pair_gap_ds_dt_other_info = -1*ones(length(sp_no)-1,2+2+1); % n-1 pairs/rows; cols: 2sp#s, distance in space, and time, on-campus
        case_no_battery = -1*ones(length(sp_no)-1,2+1); % n-1 pairs/rows; cols: 2sp#s, case#
        case_no_sleep = -1*ones(length(sp_no)-1,2+1); % n-1 pairs/rows; cols: 2sp#s, case#
        ts_bet_pair = -1*ones(length(sp_no)-1,2);
        merged_sp_seq_no = [];
        secondary_overlapped_merged_sp_pairs = [];
        secondary_overlapped_unmerged_sp_pairs = [];
        secondary_overlapped_unmerged_sp_pairs_special = []; % this is a special unmerged case where the pair can be merged as single but NOT from ref. cluster as a sequence merging
        ref_sp_centroid = [];
        for j=1:length(sp_no)-1
            ds = haversine (lat_long_sp(j,1), lat_long_sp(j,2), lat_long_sp(j+1,1), lat_long_sp(j+1,2)); % in meter bet. the median centroids
            dt = ( ts_sp(j+1,1)-ts_sp(j,2) )/(1000*60); % in minutes
            on_off = 1; % checkBothOnCampus(lat_long_sp(j:j+1,1:2)); % median centroid sps
            pair_gap_ds_dt_other_info(j,:) = [j j+1 ds dt on_off];
            ts_bet_pair(j,:) = [ts_sp(j,2) ts_sp(j+1,1)];
            %% Check and obtain Battery Charging Session Existance during data gap bet. sp-pair
            if(isempty(ts_b))
                case_no_b = 0; bet_sp_uncovered_gap_durations_b = []; new_ts_1234_b = []; %-1*ones(2,2);
            else
                [case_no_b,bet_sp_uncovered_gap_durations_b,new_ts_1234_b] = check_secondary_info_existance_bet_2sp(ts_sp(j:j+1,:),ts_b,Tmax);
            end
            case_no_battery(j,:) = [j j+1 case_no_b];
            %% Check and obtain SleepMeta Sleep Session Existance during data gap bet. sp-pair
            if(isempty(ts_s))
                case_no_s = 0; bet_sp_uncovered_gap_durations_s = []; new_ts_1234_s = []; %-1*ones(2,2);
            else
                [case_no_s,bet_sp_uncovered_gap_durations_s,new_ts_1234_s] = check_secondary_info_existance_bet_2sp(ts_sp(j:j+1,:),ts_s,Tmax);
            end
            case_no_sleep(j,:) = [j j+1 case_no_s];
            %% Apply Battery and Sleep Sessions together
            if(isempty(new_ts_1234_b) & isempty(new_ts_1234_s)) % Either NO INFO. during the gap bet. sp-pair OR NOT merged in Case-2
                if(case_no_b==2 | case_no_s==2); merged_b_s = 0; 
                else; merged_b_s = -1;
                end
                new_ts_1234_b_s = [];
            elseif(~isempty(new_ts_1234_b) & ~isempty(new_ts_1234_s)) % Both sleep and battery info. exist during gap
                dt_new_ts_1234_b = new_ts_1234_b(2,:) - new_ts_1234_b(1,:); sum_dt_b = sum(dt_new_ts_1234_b);
                dt_new_ts_1234_s = new_ts_1234_s(2,:) - new_ts_1234_s(1,:); sum_dt_s = sum(dt_new_ts_1234_s);
                dt_ext_b = new_ts_1234_b(:,2)-new_ts_1234_b(:,1);
                dt_ext_s = new_ts_1234_s(:,2)-new_ts_1234_s(:,1);
                if(sum_dt_b==0 & sum_dt_s==0) % merged in both types, i.e., battery and sleep ==> take the max. time extension
                    merged_b_s = 1;
                    if(dt_ext_b(1)>dt_ext_s(1)); new_ts_1234_b_s = new_ts_1234_b(1,:);
                    else; new_ts_1234_b_s = new_ts_1234_s(1,:);
                    end
                elseif(sum_dt_b==0) % merged in battery 
                    merged_b_s = 2;
                    new_ts_1234_b_s = new_ts_1234_b(1,:);
                elseif(sum_dt_s==0) % merged in sleep
                    merged_b_s = 3;
                    new_ts_1234_b_s = new_ts_1234_s(1,:);
                else % merged in neither battery nor sleep
                    %% take the longest extended sp#1 and sp#2 across battery and sleep
                    % sp#1
                    if(dt_ext_b(1)>dt_ext_s(1)) new_ts_1234_sp1 = new_ts_1234_b(1,:);
                    else new_ts_1234_sp1 = new_ts_1234_s(1,:);
                    end
                    % sp#2
                    if(dt_ext_b(2)>dt_ext_s(2)) new_ts_1234_sp2 = new_ts_1234_b(2,:);
                    else new_ts_1234_sp2 = new_ts_1234_s(2,:);
                    end
                    %% try to merge the longest extended sp#1 and sp#2
                    if( (new_ts_1234_sp2(1)-new_ts_1234_sp1(2))/(1000*60)<=Tmax ) % can be merged using sleep and battery together 
                        merged_b_s = 4;
                        new_ts_1234_b_s = [new_ts_1234_sp1(1) new_ts_1234_sp2(2)];
                    else
                        merged_b_s = 0;
                        new_ts_1234_b_s = [new_ts_1234_sp1; new_ts_1234_sp2];
                    end
                end
            elseif(~isempty(new_ts_1234_b) & isempty(new_ts_1234_s)) % ONLY battery info. exist during gap
                dt_new_ts_1234_b = new_ts_1234_b(2,:) - new_ts_1234_b(1,:); sum_dt_b = sum(dt_new_ts_1234_b);
                if(sum_dt_b==0) % merged in battery 
                    merged_b_s = 2;
                    new_ts_1234_b_s = new_ts_1234_b(1,:);
                else
                    merged_b_s = 0;
                    new_ts_1234_b_s = new_ts_1234_b;
                end
            elseif(isempty(new_ts_1234_b) & ~isempty(new_ts_1234_s)) % ONLY sleep info. exist during gap
                dt_new_ts_1234_s = new_ts_1234_s(2,:) - new_ts_1234_s(1,:); sum_dt_s = sum(dt_new_ts_1234_s);
                if(sum_dt_s==0) % merged in sleep
                    merged_b_s = 3;
                    new_ts_1234_b_s = new_ts_1234_s(1,:);
                else
                    merged_b_s = 0;
                    new_ts_1234_b_s = new_ts_1234_s;
                end
            end
            %% Build String of extended timestamps
            time_ext_b_s = [];
            if(merged_b_s<0) % Neither Sleep nor battery INFO available during gap
                time_ext_b_s = [];
            elseif(merged_b_s>0) % sp-pair MERGED
                time_ext_b_s = [num2str(new_ts_1234_b_s(1,1)) ';' num2str(new_ts_1234_b_s(1,2))];
            else % NOT merged, but extended
                if(~isempty(new_ts_1234_b_s))
                    time_ext_b_s = [num2str(new_ts_1234_b_s(1,1)) ';' num2str(new_ts_1234_b_s(1,2))];
                    time_ext_b_s = [time_ext_b_s ';' num2str(new_ts_1234_b_s(2,1)) ';'  num2str(new_ts_1234_b_s(2,2))];
                end
            end
            %% Build String of actual timestamps for sp#1 and sp#2
            time_act = [num2str(ts_sp(j,1)) ';' num2str(ts_sp(j,2))];
            time_act = [time_act ';' num2str(ts_sp(j+1,1)) ';'  num2str(ts_sp(j+1,2))];
            %% Write Battery Charging and SleepMeta Sleep Session's Inference on sp-pair Separately into a file 
            % Get the multiple uncovered time gaps between 2-sps for case-2 with battery info. in a string separated by ';'
            uncovered_dt_gap_b = [];
            for k=1:length(bet_sp_uncovered_gap_durations_b)
                if(k==1) uncovered_dt_gap_b = num2str(bet_sp_uncovered_gap_durations_b(k));
                else uncovered_dt_gap_b = [uncovered_dt_gap_b ';' num2str(bet_sp_uncovered_gap_durations_b(k))];
                end
            end
            % Get the timestamp extensions as strings separated by ';' when using battery charging info.
            time_ext_b = [];
            if(~isempty(new_ts_1234_b))
                time_ext_b = [num2str(new_ts_1234_b(1,1)) ';' num2str(new_ts_1234_b(1,2))];
                time_ext_b = [time_ext_b ';' num2str(new_ts_1234_b(2,1)) ';'  num2str(new_ts_1234_b(2,2))];
            end
            % Get the multiple time gaps between 2-sps for case-2 with sleep info in a string separated by ';'
            uncovered_dt_gap_s = [];
            for k=1:length(bet_sp_uncovered_gap_durations_s)
                if(k==1) uncovered_dt_gap_s = num2str(bet_sp_uncovered_gap_durations_s(k));
                else uncovered_dt_gap_s = [uncovered_dt_gap_s ';' num2str(bet_sp_uncovered_gap_durations_s(k))];
                end
            end
            % Get the timestamp extensions as strings separated by ';' when using sleep info.
            time_ext_s = [];
            if(~isempty(new_ts_1234_s))
                time_ext_s = [num2str(new_ts_1234_s(1,1)) ';' num2str(new_ts_1234_s(1,2))];
                time_ext_s = [time_ext_s ';' num2str(new_ts_1234_s(2,1)) ';'  num2str(new_ts_1234_s(2,2))];
            end
            %% Find sequence of merged sp-pairs 
            if(merged_b_s>0) 
                if(ismember(j,merged_sp_seq_no)) %% Continuation of an existing sequence
                    d_from_ref = haversine (ref_sp_centroid(1), ref_sp_centroid(2), lat_long_sp(j+1,1), lat_long_sp(j+1,2));
                    if(d_from_ref<=Dmax) %% Keep continuing sequence merging 
                        merged_sp_seq_no = [merged_sp_seq_no,j+1];
                    else %% Stop last sequence and start a new sequence 
                        ref_sp_no = 0;
                        ref_sp_centroid = [];
                        merged_sp_seq_no = [merged_sp_seq_no, -5];
                        secondary_overlapped_unmerged_sp_pairs_special = [secondary_overlapped_unmerged_sp_pairs_special; [j j+1]];
                    end
                else %% Starting a new sequence 
                    ref_sp_no = j;
                    ref_sp_centroid = lat_long_sp(j,1:2);
%                     merged_sp_seq_no = union(merged_sp_seq_no,[j j+1]);
                    merged_sp_seq_no = [merged_sp_seq_no,[j j+1]];
                    last_sp_no = j+1;
                end
            else %% pair not merged, break of a sequence 
                ref_sp_no = 0;
                ref_sp_centroid = [];
                ds_temp = haversine (lat_long_sp(j,1), lat_long_sp(j,2), lat_long_sp(j+1,1), lat_long_sp(j+1,2));
                if(merged_b_s>=0 & ds_temp<=Dmax) secondary_overlapped_unmerged_sp_pairs = [secondary_overlapped_unmerged_sp_pairs; [j j+1]]; end
            end
            % Write each of these sp-pair info into one subjLevel file
            % sp#1,sp#2,ds(meter-bet.medianCentroids),dt(minutes-bet.spBoundaries),on_off(both
            % sp on-campus?),case#(using battery),uncovered_gap_bet_sp(in minutes|using battery),case#(using sleep),
            % uncovered_gap_bet_sp(in minutes, using sleep),timeExtension(using battery),timeExtension(using sleep),
            % merged(-1=NO INFO|0=NOT Merged|1=BOTH battery and sleep merged|2=battery merged|3=sleep merged|4=battery and sleep together merged),timeExtension(2values=merged|4values=extendedButNotMerged|empty=NotExtended),actualSpPairTime
            fprintf(fileID,'%d,%d,%f,%f,%d,%d,%s,%d,%s,%s,%s,%d,%s,%s\n',j,j+1,ds,dt,on_off,case_no_b,uncovered_dt_gap_b(:),case_no_s,uncovered_dt_gap_s(:),time_ext_b(:),time_ext_s(:),merged_b_s,time_ext_b_s(:),time_act(:));
        end
        fclose(fileID);
        %% Process Sequence of Merged sp-pairs 
        if(~isempty(merged_sp_seq_no))
            merged_sp_seq_no = merged_sp_seq_no';
            merged_sp_seq_no_1st_diff = merged_sp_seq_no(2:end)-merged_sp_seq_no(1:end-1);
            merged_sp_seq_no_1st_diff = [1; merged_sp_seq_no_1st_diff];
            ind_seq_merge = find(merged_sp_seq_no_1st_diff~=1);
            seq_st_sp_no = [merged_sp_seq_no(1); merged_sp_seq_no(ind_seq_merge)];
            seq_end_sp_no = [merged_sp_seq_no(ind_seq_merge-1); merged_sp_seq_no(end)];
            secondary_overlapped_merged_sp_pairs = [seq_st_sp_no seq_end_sp_no];
            ind_v_seq = find(secondary_overlapped_merged_sp_pairs(:,1)>0);
            secondary_overlapped_merged_sp_pairs = secondary_overlapped_merged_sp_pairs(ind_v_seq,:);
        end
%         secondary_overlapped_unmerged_sp_pairs
        %% Write sequence of merged sp cluster pairs, along with unmerged sp-pairs  
        if(~isempty(secondary_overlapped_merged_sp_pairs) | ~isempty(secondary_overlapped_unmerged_sp_pairs))
            OutFile_SeqMerge = [outDir_SeqMerge 'sp_pairs_seq_pid_' num2str(pid) '.csv'];
            fileID_SeqMerge = fopen(OutFile_SeqMerge,'w+');
            fprintf(fileID_SeqMerge,'st_sp#,end_sp#,merged_unmerged(1=merged|0=unmerged||eah of these pairs has secondary info overlap during gap)\n');
            merged_unmerged = 1;
            for k=1:size(secondary_overlapped_merged_sp_pairs,1)
                fprintf(fileID_SeqMerge,'%d,%d,%d\n',secondary_overlapped_merged_sp_pairs(k,:),merged_unmerged);
            end
            merged_unmerged = 0;
            for k=1:size(secondary_overlapped_unmerged_sp_pairs,1)
                fprintf(fileID_SeqMerge,'%d,%d,%d\n',secondary_overlapped_unmerged_sp_pairs(k,:),merged_unmerged);
            end
            merged_unmerged = -1;
            for k=1:size(secondary_overlapped_unmerged_sp_pairs_special,1)
                fprintf(fileID_SeqMerge,'%d,%d,%d\n',secondary_overlapped_unmerged_sp_pairs_special(k,:),merged_unmerged);
            end
    
            fclose(fileID_SeqMerge);
        end
        
        %% General Info about the sp-pairs 
        ind_noise_affected_clustering = find(pair_gap_ds_dt_other_info(:,3)<=Dmax & pair_gap_ds_dt_other_info(:,4)<=Tmax);
        noise_affected_clustering_all_subj = [noise_affected_clustering_all_subj; [pid size(pair_gap_ds_dt_other_info,1) length(ind_noise_affected_clustering)]];
    
        %% Battery Charging Session matched sp-pairs  
        ind_BS_available_bet_sps = find(case_no_battery(:,3)>0);
        ind_BS_available_merge_sps = find(case_no_battery(:,3)==1);
        dt_merge_sps = (ts_bet_pair(ind_BS_available_merge_sps,2)-ts_bet_pair(ind_BS_available_merge_sps,1))/(1000*60); % minutes
        end_DateTime_merge_sps_sp1 = convert_timestamp_time(ts_bet_pair(ind_BS_available_merge_sps,1));
        st_DateTime_merge_sps_sp2 = convert_timestamp_time(ts_bet_pair(ind_BS_available_merge_sps,2));
        ind_midNight_sp_split = []; % this is to test percentage of sp-pairs splited just because of day-level sp-clustering approach
        for ij=1:length(ind_BS_available_merge_sps)
            temp_end = end_DateTime_merge_sps_sp1(ij,12:15);
            temp_st = st_DateTime_merge_sps_sp2(ij,12:15);
            if strcmp(temp_end,'23:5') & strcmp(temp_st,'00:0') & dt_merge_sps(ij)<=Tmax
                ind_midNight_sp_split = [ind_midNight_sp_split; ij];
            end
        end
        BS_available_all_subj = [BS_available_all_subj; [pid size(pair_gap_ds_dt_other_info,1) length(ind_BS_available_bet_sps) length(ind_BS_available_merge_sps) length(ind_midNight_sp_split)]];
        
        %% SleepMeta Sleep Session matched sp-pairs  
        ind_Slp_available_bet_sps = find(case_no_sleep(:,3)>0);
        ind_Slp_available_merge_sps = find(case_no_sleep(:,3)==1);
%         dt_merge_sps = pair_gap_ds_dt_other_info(ind_Slp_available_merge_sps,4);
        dt_merge_sps = (ts_bet_pair(ind_Slp_available_merge_sps,2)-ts_bet_pair(ind_Slp_available_merge_sps,1))/(1000*60); % minutes
        end_DateTime_merge_sps_sp1 = convert_timestamp_time(ts_bet_pair(ind_Slp_available_merge_sps,1));
        st_DateTime_merge_sps_sp2 = convert_timestamp_time(ts_bet_pair(ind_Slp_available_merge_sps,2));
%         DateTime_bet_pair = [end_DateTime_merge_sps_sp1(:,12:15) st_DateTime_merge_sps_sp2(:,12:15)]
        ind_midNight_sp_split = [];
        for ij=1:length(ind_Slp_available_merge_sps)
            temp_st = st_DateTime_merge_sps_sp2(ij,12:15);
            temp_end = end_DateTime_merge_sps_sp1(ij,12:15);
            if strcmp(temp_end,'23:5') & strcmp(temp_st,'00:0') & dt_merge_sps(ij)<=Tmax
                ind_midNight_sp_split = [ind_midNight_sp_split; ij];
            end
        end
        Slp_available_all_subj = [Slp_available_all_subj; [pid size(pair_gap_ds_dt_other_info,1) length(ind_Slp_available_bet_sps) length(ind_Slp_available_merge_sps) length(ind_midNight_sp_split)]];
    end
    return 
    
    %% Overall Noise Effect
    temp = sum(noise_affected_clustering_all_subj(:,2:3))
    temp(2)/temp(1)*100

    %% Overall Battery Charging Session matched sp-pairs
    BS_available_all_subj = double(BS_available_all_subj);
    temp = sum(BS_available_all_subj(:,2:end))
    [temp(2)/temp(1)*100 temp(3)/temp(2)*100 temp(4)/temp(2)*100]
    temp = BS_available_all_subj(:,4)./BS_available_all_subj(:,3)*100;
    figure
    boxplot(temp)
    ylabel('% merged using battery charging session')
    title('Subject Level sp-pair Merging Success')
    grid on
%     pause

    %% Overall SleepMeta Sleep Session matched sp-pairs
    Slp_available_all_subj = double(Slp_available_all_subj);
    temp = sum(Slp_available_all_subj(:,2:end))
    [temp(2)/temp(1)*100 temp(3)/temp(2)*100 temp(4)/temp(2)*100]
    temp = Slp_available_all_subj(:,4)./Slp_available_all_subj(:,3)*100;
    figure
    boxplot(temp)
    ylabel('% merged using sleep session')
    title('Subject Level sp-pair Merging Success')
    grid on
end
function [case_no,bet_sp_uncovered_gap_durations,new_ts_1234] = check_secondary_info_existance_bet_2sp(ts_sp,ts_secondary,Tmax)

    case_no = -1;
    bet_sp_uncovered_gap_durations = []; % in minutes
    
    % -1 for a st-end ts-pair means unchanged; for case-1 it'll same entries for both rows
    new_ts_1234 = -1*ones(2,2); % to capture this extended st-end time of each of the sps or merged sp; Tmax gap consideration is NEEDED !!!
    
    t1=ts_sp(1,1); t2=ts_sp(1,2);
    t3=ts_sp(2,1); t4=ts_sp(2,2);
    %% case-1
    ind1 = find(ts_secondary(:,1)<=t2 & ts_secondary(:,2)>=t3);
    if(~isempty(ind1)) % case-1 observed
%         disp('Case-1')
        case_no = 1;
        if(length(ind1)>1) % If there exist more segments by error/faulty sensor recording/processing, take the largest one
            %% Take the largest one and then compare with sp-pair timestamps to find a longer presence st-end timestamps
            dt_temp = ts_secondary(ind1,2)-ts_secondary(ind1,1); % in mili-sec
            ind_longest = find(dt_temp==max(dt_temp));
            ind1 = ind1(ind_longest);
        end
        ind1 = ind1(1); % If there exist repetitive entries, take the 1st one
        bet_sp_uncovered_gap_durations = 0;  
        %% Extensions
        % Start
        if(ts_secondary(ind1,1)<=t1) new_ts_1234(:,1) = ts_secondary(ind1,1);
        else new_ts_1234(:,1) = t1;
        end
        % End
        if(ts_secondary(ind1,2)>=t4) new_ts_1234(:,2) = ts_secondary(ind1,2);
        else new_ts_1234(:,2) = t4;
        end
    else
        %% case-2
        ind2 = find(ts_secondary(:,1)>t2 & ts_secondary(:,2)<t3);
        if(~isempty(ind2)) % case-2 observed
%             disp('Case-2')
            case_no = 2;
            %% Case-2 is quite complicated and there're not many of these. So, for these we won't go for individual cluster extensions
            % i.e., we ONLY consider the cases where sp-pair can be MERGED
            ts_sp_and_secondary_segs = [ts_sp(1,:);  ts_secondary(ind2,:); ts_sp(2,:)];
            uncovered_dt_gap_bet_2sp = ( ts_sp_and_secondary_segs(2:end,1)-ts_sp_and_secondary_segs(1:end-1,2) )./(1000*60); % in minutes
            ind_suff_cond_fail = find(uncovered_dt_gap_bet_2sp>Tmax); % SUFFICIENT condition to merge the sp-pair, i.e., each gap <= Tmax
            if(isempty(ind_suff_cond_fail))
                bet_sp_uncovered_gap_durations = uncovered_dt_gap_bet_2sp;  % in minutes
                %% Extensions -- basically the sp#1-st and sp#2-end boundaries for those merged cases ONLY
                new_ts_1234(:,1) = t1; new_ts_1234(:,2) = t4;
            end
        else % case-3 and 4 should be checked together since 2 battery seg. can overlap sp#1-ending and sp#-beginning parts separately, and the gap bet. them can be <= Tmax
            %% case-3
            ind3 = find(ts_secondary(:,2)>t2 & ts_secondary(:,2)<t3);
            %% case-4
            ind4 = find(ts_secondary(:,1)>t2 & ts_secondary(:,1)<t3);
            if(isempty(ind3) & isempty(ind4)) % NO INFO during the gap bet. sp-pair
                case_no = 0;
            elseif(~isempty(ind3) & isempty(ind4))
                case_no = 3;
                if(length(ind3)>1) 
                    %% Take the largest one and then compare with sp-pair timestamps to find a longer presence st-end timestamps
                    dt_temp = ts_secondary(ind3,2)-ts_secondary(ind3,1);
                    ind_longest = find(dt_temp==max(dt_temp));
                    ind3 = ind3(ind_longest);
                end
                ind3 = ind3(1); % If there exist repetitive entries, take the 1st one
                bet_sp_uncovered_gap_durations = (t3-ts_secondary(ind3,2))/(1000*60); % in minutes
                %% Extensions
                % Start
                if(ts_secondary(ind3,1)<=t1) new_ts_1234(:,1) = ts_secondary(ind3,1);
                else new_ts_1234(:,1) = t1;
                end
                % End
                if( (t3-ts_secondary(ind3,2))/(1000*60)<=Tmax ) new_ts_1234(:,2) = t4;
                else new_ts_1234(1,2) = ts_secondary(ind3,2); new_ts_1234(2,:) = [t3 t4];
                end
            elseif(isempty(ind3) & ~isempty(ind4))
                case_no = 4;
                if(length(ind4)>1) 
                    %% Take the largest one and then compare with sp-pair timestamps to find a longer presence st-end timestamps
                    dt_temp = ts_secondary(ind4,2)-ts_secondary(ind4,1);
                    ind_longest = find(dt_temp==max(dt_temp));
                    ind4 = ind4(ind_longest);
                end
                ind4 = ind4(1); % If there exist repetitive entries, take the 1st one
                bet_sp_uncovered_gap_durations = (ts_secondary(ind4,1)-t2)/(1000*60); % in minutes
                %% Extensions
                % End
                if(ts_secondary(ind4,2)>=t4) new_ts_1234(:,2) = ts_secondary(ind4,2);
                else new_ts_1234(:,2) = t4;
                end
                % Start
                if( (ts_secondary(ind4,1)-t2)/(1000*60)<=Tmax ) new_ts_1234(:,1) = t1;
                else new_ts_1234(1,:) = [t1 t2]; new_ts_1234(2,1) = ts_secondary(ind4,1);
                end
            elseif(~isempty(ind3) & ~isempty(ind4))
                case_no = 34;
                if(length(ind3)>1) 
                    %% Take the largest one and then compare with sp-pair timestamps to find a longer presence st-end timestamps
                    dt_temp = ts_secondary(ind3,2)-ts_secondary(ind3,1);
                    ind_longest = find(dt_temp==max(dt_temp));
                    ind3 = ind3(ind_longest);
                end
                if(length(ind4)>1) 
                    %% Take the largest one and then compare with sp-pair timestamps to find a longer presence st-end timestamps
                    dt_temp = ts_secondary(ind4,2)-ts_secondary(ind4,1);
                    ind_longest = find(dt_temp==max(dt_temp));
                    ind4 = ind4(ind_longest);
                end
                ind3 = ind3(1); % If there exist repetitive entries, take the 1st one
                ind4 = ind4(1); % If there exist repetitive entries, take the 1st one
                bet_sp_uncovered_gap_durations = (ts_secondary(ind4,1)-ts_secondary(ind3,2))/(1000*60); % in minutes
                %% Extensions
                if( (ts_secondary(ind4,1)-ts_secondary(ind3,2))/(1000*60)<=Tmax ) % merging case
                    % Start
                    if(ts_secondary(ind3,1)<=t1) new_ts_1234(:,1) = ts_secondary(ind3,1);
                    else new_ts_1234(:,1) = t1;
                    end
                    % End
                    if(ts_secondary(ind4,2)>=t4) new_ts_1234(:,2) = ts_secondary(ind4,2);
                    else new_ts_1234(:,2) = t4;
                    end
                else % NOT merging case ==> 2 separate sps with extensions
                    % 1st sp-cluster extension
                    if(ts_secondary(ind3,1)<=t1) new_ts_1234(1,:) = ts_secondary(ind3,:);
                    else new_ts_1234(1,:) = [t1 ts_secondary(ind3,2)];
                    end
                    % 2nd sp-cluster extension
                    if(ts_secondary(ind3,2)>=t4) new_ts_1234(2,:) = ts_secondary(ind4,:);
                    else new_ts_1234(2,:) = [ts_secondary(ind4,1) t4];
                    end
                end
            end
        end
    end
    %% Check if there was no timestamp extension, if so, then return empty new_ts_1234
    ind_temp = find(new_ts_1234>0);
    if(isempty(ind_temp)) new_ts_1234 = []; end
end
function str=convert_timestamp_time(curtime)
    G.TIME.TIMEZONE=-5; % Eastern Time Zone
    G.TIME.DAYLIGHTSAVING=1;
    G.TIME.FORMAT='mm/dd/yyyy HH:MM:SS';

    dnOffset = datenum('01-Jan-1970');
    tstamp=curtime/1000.0;
    timezone=G.TIME.TIMEZONE;
    dnNow = tstamp/(24*60*60) + dnOffset +timezone/24.0;
    str=datestr(dnNow,G.TIME.FORMAT);
    len=length(curtime);
    for i=1:len
            if G.TIME.DAYLIGHTSAVING==1
                if is_Daylight_Savings(str(i,:))==1
                    timezone=G.TIME.TIMEZONE+1;
                    dnNow(i) = tstamp(i)/(24*60*60) + dnOffset +timezone/24.0;
                    str(i,:)=datestr(dnNow(i),G.TIME.FORMAT);
                end
            end

    end
end
function Daylight_Savings = is_Daylight_Savings(Date)
    %Function takes a date as a string in 'mm/dd/yyyy' format and outputs a
    %logical, true if the date is during daylight savings time for that year. See
    %definition of daylight savings time in the USA
    Month = str2double(Date(1:2));
    if (Month > 3) && (Month < 11) %If the month is between April and October, true
        Daylight_Savings = true;
    elseif (Month < 3) || (Month == 12) %If month is January, February, or December, false
        Daylight_Savings = false;
    elseif Month == 3 %If month is March
        Sunday_Vect = [];
        for i = 1:31 %Loop through all 31 days of March, finding the Sundays
            if i < 10
                Day_str = sprintf('0%s',num2str(i));
            else
                Day_str = num2str(i);
            end
            March_Date = sprintf('03/%s/%s',Day_str,Date(7:10)); %String of looping march date
            DOW = weekday(March_Date,'mm/dd/yyyy'); 
            if DOW == 1 %DOW is 1 if the day is a Sunday
                Sunday_Vect = [Sunday_Vect i];
            end
        end
        DS_Day = Sunday_Vect(2); %Take the 2nd Sunday in March
        My_Day = str2double(Date(4:5));
        if My_Day >= DS_Day %If on or after 2nd Sunday, true
            Daylight_Savings = true; 
        else
            Daylight_Savings = false;
        end
    elseif Month == 11 %If month is November
        Sunday_Vect = [];
        for i = 1:30 %Loop through all 30 days of November, finding the Sundays
            if i < 10
                 Day_str = sprintf('0%s',num2str(i));
            else
                Day_str = num2str(i);
            end
            Nov_Date = sprintf('11/%s/%s',Day_str,Date(7:10)); %String of looping november date
            DOW = weekday(Nov_Date,'mm/dd/yyyy'); 
            if DOW == 1 %DOW is 1 if the day is a Sunday
                Sunday_Vect = [Sunday_Vect i];
            end
        end
        DS_Day = Sunday_Vect(1); %Take the 1st Sunday of the month
        My_Day = str2double(Date(4:5));
        if My_Day < DS_Day %If it is before the 1st Sunday, true
            Daylight_Savings = true;
        else
            Daylight_Savings = false;
        end
    end
end