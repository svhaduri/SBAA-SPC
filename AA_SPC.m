function AA_SPC()

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
    
    disp('---------------------------------------------------------------------')
    disp('Please enter ')
    disp('1 to execute evaluating_activity_assisted_spc_netHealth()')
    disp('anything else to execute validating_activity_assisted_spc_netHealth()')
    disp('---------------------------------------------------------------------')
    op=input('Your choice : ');
    
    %% USEFUL OUTDOOR ('0', i.e., 'NO POI') walk Sessions' Activity minutes
    File = [Dir 'inDir/' 'useful_oncampus_outdoor_Fitbit_LocAvailable_labelled_walk_sessions' '.csv'];
    [pid_wlk_out,ts_wlk_out,total_steps_wlk_out]=read_useful_wlk_sessions(File);
    %% USEFUL INDOOR (i.e., 'non-zero POI') walk Sessions' Activity minutes 
    File = [Dir 'inDir/' 'useful_oncampus_indoor_Fitbit_LocAvailable_labelled_walk_sessions' '.csv'];
    [pid_wlk_in,ts_wlk_in,total_steps_wlk_in]=read_useful_wlk_sessions(File);
    %% Threshold Calculation 
    dt_out = (ts_wlk_out(:,2)-ts_wlk_out(:,1))/(1000*60); dt_in = (ts_wlk_in(:,2)-ts_wlk_in(:,1))/(1000*60); % in minutes
    avg_steps_val = [total_steps_wlk_out./dt_out;  total_steps_wlk_in./dt_in]; avg_steps_grp = [ones(length(total_steps_wlk_out),1); 2*ones(length(total_steps_wlk_in),1)];
    disp('Statistics of Avg. Step Count (per minute) during walk sessions ....')
    [mean(avg_steps_val) std(avg_steps_val) std(avg_steps_val)/sqrt(length(avg_steps_val))]
    outlier_range_wlk_out_in = [prctile(avg_steps_val,25)-1.5*iqr(avg_steps_val) prctile(avg_steps_val,75)+1.5*iqr(avg_steps_val)]
    pid_wlk = [pid_wlk_out; pid_wlk_in];
    ts_wlk = [ts_wlk_out; ts_wlk_in];
    dt_wlk = [dt_out; dt_in];
    
    ind_v_subj = find(ismember(pid_wlk,subjList));
    pid_wlk = pid_wlk(ind_v_subj,:);
    ts_wlk = ts_wlk(ind_v_subj,:);
%     unique(pid_wlk)
%     pause
    
    if(op==1) %% Activity-Assisted Stay Point Clustering on UNLABELLED Dataset
        evaluating_activity_assisted_spc_netHealth(pid_wlk,ts_wlk); 
    else %% Activity-Assisted Stay Point Clustering (AA-SPC) on Ground Truth (i.e., LABELLED) Dataset
        validating_activity_assisted_spc_netHealth(pid_wlk,ts_wlk); % w = 10 minutes = Tmax
    end
end
function validating_activity_assisted_spc_netHealth(pid_wlk,ts_wlk)
    global Dir    
    global subjList
    global Tmin Tmax Dmax
    
    POP_POA_POU_sp_pair_count_all_subj=[]; % Proof-Of-Presence, Proof-Of-Absence, and Proof-Of-Undecided 
    POP_POA_cluster_level_stay_duration=[]; POP_cluster_level_stay_duration=[];
    POP_POA_subj_level_avg_stay_duration_diff=[]; POP_subj_level_avg_stay_duration_diff=[];
    POP_total_sp_duration=[];
    candidate_pids_sp_pair_D_C=[];
    step_counts_wlk_all_subj=[];
    parfor i=1:length(subjList)
%     for i=1:length(subjList)
        pid = subjList(i); 
        pid
        
        %% Read the USEFUL REPRESENTATIVE sp-clusters (>=60 minutes long, activity data available during 2 20-minute endings and in between them) 
        File = [Dir 'inDir/SubjectLevel_REPRESENTATIVE_minDuration60minutes_sp_median_Tmin' num2str(Tmin) 'minutes_Dmax' num2str(Dmax) 'meters_with_Place_Labels/' 'representative_useful_sp_pid_' num2str(pid) '.csv'];
        if(~exist(File,'file')); disp([File ' Does''t exist!!']); continue; end
        fid2 = fopen(File);
        % sp#,st_ts,end_ts,BldgID
        B = textscan(fid2,'%d%f%f%d','delimiter',',','headerlines',1);
        fclose(fid2);
        cno_c=B{1,1}; % 'cno_c' cno = cluster no, 'c' subscript = cluster
        if(length(cno_c)<2); continue; end % at least 2 sps are needed to form pair
        [st_ts_c,end_ts_c]=B{1,2:3}; % mili-sec
        BldgID=B{1,4}; % majority voted single building IDs that are matched across rectangle and radius approaches
        
        %% Find sp pairs from the same bldg
        poiTypeIDs=[1:31, 36:74, 76:86, 95:98]; 
        act_ind_succ_sp=find(ismember(BldgID,poiTypeIDs'));
        if(isempty(act_ind_succ_sp)); disp([num2str(pid) ' doesn''t have a sp from the bldg that will be splited later for synthetic data-based validation of AA-SPC ...']); continue; end
        candidate_pids_sp_pair_D_C=[candidate_pids_sp_pair_D_C;pid];
%         pause

        %% Read raw valid activity data (i.e., heart rate available)
        File = [Dir 'inDir/SubjectLevel_Activity_Minutes_Merged_ValidMinutes/' 'NetHealth_' num2str(pid) '_Activity' '.csv'];
        if(~exist(File,'file')); continue; end
        fid2 = fopen(File);
        % ts(mili-sec),fitStepCount
        B = textscan(fid2,'%f%d','delimiter',',','headerlines',1);
        fclose(fid2);
        [ts_f_a,steps_f_a]=B{1,:}; % subscript 'f' and 'a' stand for fitbit and activity
        
        %% Compute mean and SD of step counts during walk sessions of this particular pid
        ind_wlk=find(pid_wlk==pid);
        step_counts_wlk=[];
        for j=1:length(ind_wlk)
            ind_acti=find(ts_f_a>=ts_wlk(ind_wlk(j),1) & ts_f_a<=ts_wlk(ind_wlk(j),2));
            if(isempty(ind_acti)); continue; end
            step_counts_wlk=[step_counts_wlk,double(steps_f_a(ind_acti))'];
        end
        if(isempty(step_counts_wlk)); disp([num2str(pid) ' doesn''t have activity sample during walk sessions >>>']); continue; end
        step_counts_wlk_all_subj=[step_counts_wlk_all_subj,step_counts_wlk];
        mean_wlk=mean(step_counts_wlk); std_wlk=std(step_counts_wlk); 
%         [mean_wlk std_wlk]
        range_wlk_personalized=[mean_wlk-std_wlk mean_wlk+std_wlk];
        
        %% Peronalized wlk band/range
        range_wlk=range_wlk_personalized;
        
        %% Check whether sp pairs from the same bldg can be merged using activity (i.e.,  step count) information
        POP_sp_pair_ind=[]; % Proof-Of-Presence, i.e., mergable sp-pairs indices
        POA_sp_pair_ind=[]; % Proof-Of-Absence, i.e., NOT mergable sp-pairs indices
        POU_sp_pair_ind=[]; % Proof-Of-Undecided sp-pairs indices
        POP_sp_duration=[];
        POA_sp_duration=[]; 
%         size(act_ind_succ_sp,1)
%         pause
        for j=1:size(act_ind_succ_sp,1)
%             j
            t1=st_ts_c(act_ind_succ_sp(j,1));  t2=t1+(20*60*1000); % 20 minute at the begining to mimic as sp1
            t4=end_ts_c(act_ind_succ_sp(j,1)); t3=t4-(20*60*1000); % 20 minute at the end to mimic as sp2 
            
            %% Compute mean and SD of step counts during s1 and s2 together for this particular sp-pair of this pid
            ind_acti_s1=find(ts_f_a>=t1 & ts_f_a<=t2); ind_acti_s2=find(ts_f_a>=t3 & ts_f_a<=t4); ind_acti_s1s2=[ind_acti_s1;ind_acti_s2];
            mean_s1s2=mean(double(steps_f_a(ind_acti_s1s2))); std_s1s2=std(double(steps_f_a(ind_acti_s1s2))); 
            range_s1s2=[mean_s1s2-std_s1s2 mean_s1s2+std_s1s2];
        
            t=t2;
            w=10*60*1000; % 10 minutes window
            flag_POP=1; % Proof-Of-Presence(POP)=1
            while t<t3
                w_st=t; w_end=t+w;
                if(w_end>t3) % windows shorter than w = 10 minutes are fine to tolerate since Tmax = 10 minutes
                    w_end=t3; 
                    break;
                end
                ind_acti_win=find(ts_f_a>=w_st & ts_f_a<=w_end);
                %% First, test activity values with respect to walk range
                if(isempty(ind_acti_win)); flag_POP=-2; break; end % UNDEFINED -- NO activity (i.e., step count) data between for a complete window
                ind_wlk_range=find(steps_f_a(ind_acti_win)>=range_wlk(1) & steps_f_a(ind_acti_win)<=range_wlk(2));
                if(length(ind_wlk_range)/length(ind_acti_win)>=0.5) % Majority voting to decide POA
                    POA_sp_pair_ind=[POA_sp_pair_ind; act_ind_succ_sp(j,:)];
                    POA_sp_duration=[POP_sp_duration; [40*60*1000 40*60*1000]];
                    flag_POP=0;
                    break; 
                end
                %% Now, test activity values with respect to s1 and s2's activity range
                if(isempty(ind_acti_s1s2)); break; end % UNDEFINED -- NO activity (i.e., step count) data during s1 and s2
                ind_s1s2_range=find(steps_f_a(ind_acti_win)>=range_s1s2(1) & steps_f_a(ind_acti_win)<=range_s1s2(2));
                if(length(ind_s1s2_range)/length(ind_acti_win)<0.5)  % Majority voting to decide POU
                    POU_sp_pair_ind=[POU_sp_pair_ind; act_ind_succ_sp(j,:)];
                    flag_POP=-1;
                    break; 
                else % So far so good ==> update window time counter
                    t=w_end;
                end
            end
            if(flag_POP==1) 
                POP_sp_pair_ind=[POP_sp_pair_ind; act_ind_succ_sp(j,:)]; 
                POP_sp_duration=[POP_sp_duration; [40*60*1000 (t4-t1)]];
                POP_total_sp_duration=[POP_total_sp_duration; [(t4-t1)/(1000*60) (t3-t2)/(1000*60)]]; % in minutes
            end
        end
%         POP_sp_pair_ind
%         POA_sp_pair_ind
%         POU_sp_pair_ind
        POP_POA_cluster_level_stay_duration=[];
    
        POP_POA_POU_sp_pair_count_all_subj=[POP_POA_POU_sp_pair_count_all_subj; [pid size(POP_sp_pair_ind,1) size(POA_sp_pair_ind,1) size(POU_sp_pair_ind,1)]];
        %% Duration calculation considering both POP and POA
        POP_POA_cluster_level_stay_duration=[POP_POA_cluster_level_stay_duration; POP_sp_duration; POA_sp_duration]; 
        %% Duration calculation considering ONLY POP
        POP_cluster_level_stay_duration=[POP_cluster_level_stay_duration; POP_sp_duration]; 
        percent_POP_cluster_level_stay_duration_increase=[]; percent_POA_cluster_level_stay_duration_increase=[];
        if(~isempty(POP_sp_duration)); percent_POP_cluster_level_stay_duration_increase=(POP_sp_duration(:,2)-POP_sp_duration(:,1))./POP_sp_duration(:,2)*100; end
        if(~isempty(POA_sp_duration)); percent_POA_cluster_level_stay_duration_increase=(POA_sp_duration(:,2)-POA_sp_duration(:,1))./POA_sp_duration(:,2)*100; end
        %% Duration calculation considering both POP and POA
        POP_POA_subj_level_avg_stay_duration_diff=[POP_POA_subj_level_avg_stay_duration_diff; [pid mean([percent_POP_cluster_level_stay_duration_increase; percent_POA_cluster_level_stay_duration_increase])]];
        %% Duration calculation considering ONLY POP
        POP_subj_level_avg_stay_duration_diff=[POP_subj_level_avg_stay_duration_diff; [pid mean(percent_POP_cluster_level_stay_duration_increase)]];
    
    end
%     POP_POA_POU_sp_pair_count_all_subj
    size(POP_POA_POU_sp_pair_count_all_subj) % 368 subjects with USEFUL REPRESENTATIVE sp-clusters
    sum(POP_POA_POU_sp_pair_count_all_subj)
%     pause
    %% Analysis of % clusters successfully merged using AA-SPC 
    POP_POA_POU_sp_pair_count_all_subj=double(POP_POA_POU_sp_pair_count_all_subj);
    %% without considering POU
%     percent_merged=POP_POA_POU_sp_pair_count_all_subj(:,2)./(POP_POA_POU_sp_pair_count_all_subj(:,2)+POP_POA_POU_sp_pair_count_all_subj(:,3))*100; 
    %% considering POU in addition to POA and POP
    percent_merged=POP_POA_POU_sp_pair_count_all_subj(:,2)./sum(POP_POA_POU_sp_pair_count_all_subj(:,2:4),2)*100; 
    figure
    subplot(1,2,1)
    disp('Mean SD of % merged ...')
    [mean(percent_merged) std(percent_merged)]
    boxplot(percent_merged,'Labels',{'(a)'})
    ylabel('% merged')
    grid on
    %% Analysis of % clusters count decreased using AA-SPC for subject-level
    %% without considering POU
%     CC_before_merging=2*sum(POP_POA_POU_sp_pair_count_all_subj(:,2:3),2); 
%     CC_after_merging=POP_POA_POU_sp_pair_count_all_subj(:,2)+2*POP_POA_POU_sp_pair_count_all_subj(:,3);
    %% considering POU in addition to POA and POP
    CC_before_merging=2*sum(POP_POA_POU_sp_pair_count_all_subj(:,2:4),2); 
    CC_after_merging=POP_POA_POU_sp_pair_count_all_subj(:,2)+2*sum(POP_POA_POU_sp_pair_count_all_subj(:,3:4),2);
    CCD=(CC_before_merging-CC_after_merging)./CC_before_merging*100;
%     figure
    subplot(1,2,2)
    disp('Mean SD of CCD ...')
    [mean(CCD) std(CCD)]
    boxplot(CCD,'Labels',{'(b)'})
    ylabel('% CCD')
    grid on
    %% Analysis of % clusters times increased using AA-SPC both for cluster-level across all subjects together and subject avg.
    figure
    subplot(1,2,1)
    %% Duration calculation considering both POP and POA
%     percent_POP_POA_cluster_level_stay_duration_increase=(POP_POA_cluster_level_stay_duration(:,2)-POP_POA_cluster_level_stay_duration(:,1))./POP_POA_cluster_level_stay_duration(:,2)*100;
    %% Duration calculation considering ONLY POP
    percent_POP_POA_cluster_level_stay_duration_increase=(POP_cluster_level_stay_duration(:,2)-POP_cluster_level_stay_duration(:,1))./POP_cluster_level_stay_duration(:,2)*100;
    disp('Paired T-test on % CTI ...')
    [h,p,ci,stats] = ttest(percent_POP_POA_cluster_level_stay_duration_increase)
    disp('Cluster-level Mean SD of CTI ...')
    [mean(percent_POP_POA_cluster_level_stay_duration_increase) std(percent_POP_POA_cluster_level_stay_duration_increase)]
    boxplot(percent_POP_POA_cluster_level_stay_duration_increase,'Labels',{'(a)'})
    ylabel('% CTI')
    grid on
%     figure
    subplot(1,2,2)
    %% Duration calculation considering both POP and POA
%     POP_POA_subj_level_avg_stay_duration_diff=double(POP_POA_subj_level_avg_stay_duration_diff);
    %% Duration calculation considering ONLY POP
    POP_POA_subj_level_avg_stay_duration_diff=double(POP_subj_level_avg_stay_duration_diff);
    disp('Subject-level Mean SD of CTI ...')
    [mean(POP_POA_subj_level_avg_stay_duration_diff(:,2)) std(POP_POA_subj_level_avg_stay_duration_diff(:,2))]
    boxplot(POP_POA_subj_level_avg_stay_duration_diff(:,2),'Labels',{'(b)'})
    ylabel('Avg. CTI (%)')
    grid on
%     pause
    %% PDF and CDF of CTI cluster-level across all subjects
    inc = 10; %5 %  
    dt_improve_count_distribution = zeros(100/inc,1);
    bin_counter = 0;
    for j=0:inc:100-inc
        bin_counter = bin_counter + 1;
        ind_bin_count = find(percent_POP_POA_cluster_level_stay_duration_increase>j & percent_POP_POA_cluster_level_stay_duration_increase<=j+inc);
        dt_improve_count_distribution(bin_counter,1) = length(ind_bin_count);
    end
    p = dt_improve_count_distribution/sum(dt_improve_count_distribution);
    c = cumsum(p);
    figure
    b = bar(p,'c');
    xlim([.5 100/inc+.5])
    xtickLabels = {'10','20','30','40','50','60','70','80','90','100'};
    xticklabels(xtickLabels)
    hold on 
    plot(c,'b--'); 
    xlabel('% Increase in Cluster Time')
    ylabel('Probability')
    title('AA-SPC vs. ESPC')
    grid on
    %% Analysis of cluster time/duration and data gap of successfully merged clusters
    disp('Avg and SD of Duration (in minutes) of successfulluy merged clusters')
    [mean(POP_total_sp_duration(:,1)) std(POP_total_sp_duration(:,1))]
    figure 
    boxplot(POP_total_sp_duration(:,1))
    grid on
    ylabel('Cluster Duration (minuites)')
    disp('Avg and SD of data gap (in minutes) between pairs of successfulluy merged clusters')
    [mean(POP_total_sp_duration(:,2)) std(POP_total_sp_duration(:,2))]
    disp('Min Max of data gap (in minutes) between pairs of successfulluy merged clusters')
    [min(POP_total_sp_duration(:,2)) max(POP_total_sp_duration(:,2))]
    figure 
    boxplot(POP_total_sp_duration(:,2))
    grid on
    ylabel('Length of Data Gap (minutes)')
    %% PDF and CDF of data gaps of successfully merged clusters 
    inc = 20;   
    upper_limit_gap_duration=1400; %max(POP_total_sp_duration(:,2)); % 1400 minutes = 24 hour
    data_gap_count_distribution = zeros(upper_limit_gap_duration/inc,1);
%     data_gap_count_distribution = zeros(1400/inc,1);
    bin_counter = 0;
    for j=0:inc:upper_limit_gap_duration-inc
        bin_counter = bin_counter + 1;
        ind_bin_count = find(POP_total_sp_duration(:,2)>j & POP_total_sp_duration(:,2)<=j+inc);
        data_gap_count_distribution(bin_counter,1) = length(ind_bin_count);
    end
    disp(['Total Number of Successfully Merged Cluster Pairs considered for PDF/CDF : ' num2str(sum(data_gap_count_distribution))])
    p = data_gap_count_distribution/sum(data_gap_count_distribution);
    c = cumsum(p);
    ind_90=find(c>=.9);
    figure
    b = bar(p,'c');
    ylim([0 1])
    xlim([.5 upper_limit_gap_duration/inc+.5])
    hold on 
    plot(c,'bo--'); 
    plot([0 ind_90(1)],[c(ind_90(1)) c(ind_90(1))],'k-'); hold on; % horizontal line at cdf = .9
    plot([ind_90(1) ind_90(1)],[0 c(ind_90(1))],'k:'); hold on; % vertical line at cdf = .9
    xlabel('Length of Data Gap (minutes)')
    ylabel('Probability')
    grid on
end
function evaluating_activity_assisted_spc_netHealth(pid_wlk,ts_wlk)
    global Dir
    global subjList
    global Tmin Tmax Dmax

    outDir_SeqMerge = [Dir 'outDir/SubjectLevel_sp_median_pairs_Tmin' num2str(Tmin) 'minutes_Dmax' num2str(Dmax) 'meters_SequenceMerging_AASPC/'];
    if(~exist(outDir_SeqMerge,'dir')); mkdir(outDir_SeqMerge); end
    
    POP_POA_POU_sp_pair_count_all_subj=[]; % Proof-Of-Presence, Proof-Of-Absence, and Proof-Of-Undecided 
    POP_POA_sp_cluster_count_all_subj=[]; POP_POA_POU_sp_cluster_count_all_subj=[];
    POP_POA_subject_level_stay_duration_all_subj=[];
    POP_POA_cluster_level_stay_duration=[]; POP_POA_subj_level_avg_stay_duration_diff=[];
    candidate_pids_sp_pair_D_C=[];
    step_counts_wlk_all_subj=[];
    parfor i=1:length(subjList)
%     for i=1:length(subjList)
        pid = subjList(i); 
        pid
        
        %% Read the USEFUL (i.e., same label from rectangle and radius approaches) sp-clusters 
        File = [Dir 'inDir/SubjectLevel_USEFUL_sp_median_Tmin' num2str(Tmin) 'minutes_Dmax' num2str(Dmax) 'meters_with_Place_Labels/' 'useful_sp_pid_' num2str(pid) '.csv'];
        if(~exist(File,'file')); disp([File ' Does''t exist!!']); continue; end
        fid2 = fopen(File);
        % sp#,st_ts,end_ts,BldgID
        B = textscan(fid2,'%d%f%f%d','delimiter',',','headerlines',1);
        fclose(fid2);
        cno_c=B{1,1}; % 'cno_c' cno = cluster no, 'c' subscript = cluster
        if(length(cno_c)<2); continue; end % at least 2 sps are needed to form pair
        [st_ts_c,end_ts_c]=B{1,2:3}; % mili-sec
        ts_c = [st_ts_c end_ts_c];
        BldgID=B{1,4}; % majority voted single building IDs that are matched across rectangle and radius approaches
        
        %% Find sp pairs from the same bldg
        act_ind_succ_sp=[];
        poiTypeIDs=[1:31,36:74, 76:86, 95:98]; 
        for j=1:length(poiTypeIDs)
            bldg_id=poiTypeIDs(j);
            ind_D_C_c=find(BldgID==bldg_id);
            if(isempty(ind_D_C_c)); continue; end
            diff_sp_no=cno_c(ind_D_C_c(2:end))-cno_c(ind_D_C_c(1:end-1)); 
            temp_ind_succ_sp=find(diff_sp_no==1);
            for k=1:length(temp_ind_succ_sp)
                act_ind_succ_sp=[act_ind_succ_sp;[ind_D_C_c(temp_ind_succ_sp(k)) ind_D_C_c(temp_ind_succ_sp(k)+1) bldg_id]];
            end
        end
%         act_ind_succ_sp
        if(isempty(act_ind_succ_sp)); disp([num2str(pid) ' doesn''t have a pair of sps from the same bldg ...']); continue; end
        candidate_pids_sp_pair_D_C=[candidate_pids_sp_pair_D_C;pid];
%         pause

        %% Read raw valid activity data (i.e., heart rate available)
        File = [Dir 'inDir/SubjectLevel_Activity_Minutes_Merged_ValidMinutes/' 'NetHealth_' num2str(pid) '_Activity' '.csv'];
        if(~exist(File,'file')); continue; end
        fid2 = fopen(File);
        % ts(mili-sec),fitStepCount
        B = textscan(fid2,'%f%d','delimiter',',','headerlines',1);
        fclose(fid2);
        [ts_f_a,steps_f_a]=B{1,:}; % subscript 'f' and 'a' stand for fitbit and activity
        
        %% Compute mean and SD of step counts during walk sessions of this particular pid
        ind_wlk=find(pid_wlk==pid);
        step_counts_wlk=[];
        for j=1:length(ind_wlk)
            ind_acti=find(ts_f_a>=ts_wlk(ind_wlk(j),1) & ts_f_a<=ts_wlk(ind_wlk(j),2));
            if(isempty(ind_acti)); continue; end
            step_counts_wlk=[step_counts_wlk,double(steps_f_a(ind_acti))'];
        end
        if(isempty(step_counts_wlk)); disp([num2str(pid) ' doesn''t have activity sample during walk sessions >>>']); continue; end
        step_counts_wlk_all_subj=[step_counts_wlk_all_subj,step_counts_wlk];
        mean_wlk=mean(step_counts_wlk); std_wlk=std(step_counts_wlk); 
%         [mean_wlk std_wlk]
        range_wlk_personalized=[mean_wlk-std_wlk mean_wlk+std_wlk];
        
        %% Peronalized wlk band/range
        range_wlk=range_wlk_personalized;
        
        %% Check whether sp pairs from the same bldg can be merged using activity (i.e.,  step count) information
        POP_sp_pair_ind=[]; % Proof-Of-Presence, i.e., mergable sp-pairs indices
        POA_sp_pair_ind=[]; % Proof-Of-Absence, i.e., NOT mergable sp-pairs indices
        POU_sp_pair_ind=[]; % Proof-Of-Undecided sp-pairs indices
        for j=1:size(act_ind_succ_sp,1)
%             j
            t1=st_ts_c(act_ind_succ_sp(j,1)); t2=end_ts_c(act_ind_succ_sp(j,1));
            t3=st_ts_c(act_ind_succ_sp(j,2)); t4=end_ts_c(act_ind_succ_sp(j,2));
            
            %% Compute mean and SD of step counts during s1 and s2 together for this particular sp-pair of this pid
            ind_acti_s1=find(ts_f_a>=t1 & ts_f_a<=t2); ind_acti_s2=find(ts_f_a>=t3 & ts_f_a<=t4); ind_acti_s1s2=[ind_acti_s1;ind_acti_s2];
            mean_s1s2=mean(double(steps_f_a(ind_acti_s1s2))); std_s1s2=std(double(steps_f_a(ind_acti_s1s2))); 
            range_s1s2=[mean_s1s2-std_s1s2 mean_s1s2+std_s1s2];
        
            t=t2;
            w=10*60*1000; % 10 minutes window
            flag_POP=1; % Proof-Of-Presence(POP)=1
            while t<t3
                w_st=t; w_end=t+w;
                if(w_end>t3) % windows shorter than w = 10 minutes are fine to tolerate since Tmax = 10 minutes
                    w_end=t3; 
                    break;
                end
                ind_acti_win=find(ts_f_a>=w_st & ts_f_a<=w_end);
                %% First, test activity values with respect to walk range
                if(isempty(ind_acti_win)); break; end % UNDEFINED -- NO activity (i.e., step count) data between for a complete window
                ind_wlk_range=find(steps_f_a(ind_acti_win)>=range_wlk(1) & steps_f_a(ind_acti_win)<=range_wlk(2));
                if(length(ind_wlk_range)/length(ind_acti_win)>=0.5) % Majority voting to decide POA
                    POA_sp_pair_ind=[POA_sp_pair_ind; act_ind_succ_sp(j,:)];
                    flag_POP=0;
                    break; 
                end
                %% Now, test activity values with respect to s1 and s2's activity range
                if(isempty(ind_acti_s1s2)); flag_POP=-2; break; end % UNDEFINED -- NO activity (i.e., step count) data during s1 and s2
                ind_s1s2_range=find(steps_f_a(ind_acti_win)>=range_s1s2(1) & steps_f_a(ind_acti_win)<=range_s1s2(2));
                if(length(ind_s1s2_range)/length(ind_acti_win)<0.5)  % Majority voting to decide POU
                    POU_sp_pair_ind=[POU_sp_pair_ind; act_ind_succ_sp(j,:)];
                    flag_POP=-1;
                    break; 
                else % So far so good ==> update window time counter
                    t=w_end;
                end
            end
            if(flag_POP==1); POP_sp_pair_ind=[POP_sp_pair_ind; act_ind_succ_sp(j,:)]; end
        end
%         POP_sp_pair_ind
%         POA_sp_pair_ind
%         POU_sp_pair_ind
        POP_POA_POU_sp_pair_count_all_subj=[POP_POA_POU_sp_pair_count_all_subj; [pid size(POP_sp_pair_ind,1) size(POA_sp_pair_ind,1) size(POU_sp_pair_ind,1)]];
        
        %% Sequence merging for POP sp-pairs 
        POP_sp_pair_ind_seq=[];
%         size(POP_sp_pair_ind)
        u_D_clss_bldg_ids=unique(POP_sp_pair_ind(:,end));
        for j=1:length(u_D_clss_bldg_ids)
%             u_D_clss_bldg_ids(j)
            temp_pair_ind=find(POP_sp_pair_ind(:,end)==u_D_clss_bldg_ids(j));
            seq_st_ind=POP_sp_pair_ind(temp_pair_ind(1),1);
            for k=1:length(temp_pair_ind)-1
%                 k
                if(POP_sp_pair_ind(temp_pair_ind(k),2)~=POP_sp_pair_ind(temp_pair_ind(k+1),1))
                    POP_sp_pair_ind_seq=[POP_sp_pair_ind_seq;[seq_st_ind POP_sp_pair_ind(temp_pair_ind(k),2) u_D_clss_bldg_ids(j)]];
                    seq_st_ind=POP_sp_pair_ind(temp_pair_ind(k+1),1);
                end
            end
            POP_sp_pair_ind_seq=[POP_sp_pair_ind_seq;[seq_st_ind POP_sp_pair_ind(temp_pair_ind(end),2) u_D_clss_bldg_ids(j)]];
%             pause
        end

        %% Write sequence of merged sp cluster pairs, along with unmerged sp-pairs  
        [s,IX]=sort(POP_sp_pair_ind_seq(:,1));
        POP_sp_pair_ind_seq=POP_sp_pair_ind_seq(IX,:);
        if(~isempty(POP_sp_pair_ind_seq))
            OutFile_SeqMerge = [outDir_SeqMerge 'sp_pairs_seq_pid_' num2str(pid) '.csv'];
            fileID_SeqMerge = fopen(OutFile_SeqMerge,'w+');
            fprintf(fileID_SeqMerge,'st_sp#,end_sp#,merged_unmerged(1=merged|0=unmerged|-1=undecided||each of these pairs has step counts during gap)\n');
            merged_unmerged = 1; % Merged cases after sequence merging, i.e., POP
            for k=1:size(POP_sp_pair_ind_seq,1)
                fprintf(fileID_SeqMerge,'%d,%d,%d\n',cno_c(POP_sp_pair_ind_seq(k,1)),cno_c(POP_sp_pair_ind_seq(k,2)),merged_unmerged);
            end
            merged_unmerged = 0; % Unmerged cases, i.e., POA
            for k=1:size(POA_sp_pair_ind,1)
                fprintf(fileID_SeqMerge,'%d,%d,%d\n',cno_c(POA_sp_pair_ind(k,1)),cno_c(POA_sp_pair_ind(k,2)),merged_unmerged);
            end
            merged_unmerged = -1; % Undecided cases, i.e., POU
            for k=1:size(POU_sp_pair_ind,1)
                fprintf(fileID_SeqMerge,'%d,%d,%d\n',cno_c(POU_sp_pair_ind(k,1)),cno_c(POU_sp_pair_ind(k,2)),merged_unmerged);
            end
            fclose(fileID_SeqMerge);
        end 
%         continue;
%         pause
        
        %% Cluster Count before-after AA-SPC
        single_POP_ind=[]; single_POA_ind=[]; single_POU_ind=[]; 
        if(~isempty(POP_sp_pair_ind)); single_POP_ind=union(POP_sp_pair_ind(:,1),POP_sp_pair_ind(:,2)); end
        if(~isempty(POA_sp_pair_ind)); single_POA_ind=union(POA_sp_pair_ind(:,1),POA_sp_pair_ind(:,2)); end
        if(~isempty(POU_sp_pair_ind)); single_POU_ind=union(POU_sp_pair_ind(:,1),POU_sp_pair_ind(:,2)); end
        single_POP_POA_ind=union(single_POP_ind,single_POA_ind);
        single_POP_POA_POU_ind=union(single_POP_POA_ind,single_POU_ind);
        count_POA_POP_before_merging=length(single_POP_POA_ind);
        count_POA_POP_POU_before_merging=length(single_POP_POA_POU_ind);
        single_POP_ind_after_seq_mergeing=[];
        for j=1:size(POP_sp_pair_ind_seq,1)
            single_POP_ind_after_seq_mergeing=[single_POP_ind_after_seq_mergeing, POP_sp_pair_ind_seq(j,1):POP_sp_pair_ind_seq(j,2)];
        end
        single_POA_ind_after_seq_merged=setdiff(single_POP_POA_ind,single_POP_ind_after_seq_mergeing');
        count_POP_after_seq_merging=size(POP_sp_pair_ind_seq,1);
        count_POA_after_seq_merging=length(single_POA_ind_after_seq_merged);
        count_POA_POP_after_merging=count_POP_after_seq_merging+count_POA_after_seq_merging;
        POP_POA_sp_cluster_count_all_subj=[POP_POA_sp_cluster_count_all_subj;[pid count_POA_POP_before_merging count_POA_POP_after_merging]]; % previous calculation of CCD didn't include POU
        single_POA_POU_ind_after_seq_merged=setdiff(single_POP_POA_POU_ind,single_POP_ind_after_seq_mergeing');
        count_POA_POU_after_seq_merging=length(single_POA_POU_ind_after_seq_merged);
        count_POA_POP_POU_after_merging=count_POP_after_seq_merging+count_POA_POU_after_seq_merging;
        POP_POA_POU_sp_cluster_count_all_subj=[POP_POA_POU_sp_cluster_count_all_subj;[pid count_POA_POP_POU_before_merging count_POA_POP_POU_after_merging]]; % new way of calculating CCD is based on POP, POA, and POU
        
        %% Cluster Time before-after AA-SPC on a subject-level
        stay_duration_before_merging=[];
        for j=1:length(single_POP_ind) % ONLY POP 
%         for j=1:length(single_POP_POA_ind) % POP and POA together 
            stay_duration_before_merging=[stay_duration_before_merging;(ts_c(single_POP_POA_ind(j),2)-ts_c(single_POP_POA_ind(j),1))];
        end
        stay_duration_after_merging=[];
        for j=1:size(POP_sp_pair_ind_seq,1) % Sequence mereged clusters 
            stay_duration_after_merging=[stay_duration_after_merging;(ts_c(POP_sp_pair_ind_seq(j,2),2)-ts_c(POP_sp_pair_ind_seq(j,1),1))];
        end
%         for j=1:length(single_POA_ind_after_seq_merged) % Unmerged clusters
%             stay_duration_after_merging=[stay_duration_after_merging;(ts_c(single_POA_ind_after_seq_merged(j),2)-ts_c(single_POA_ind_after_seq_merged(j),1))];
%         end
        POP_POA_subject_level_stay_duration_all_subj=[POP_POA_subject_level_stay_duration_all_subj;[pid sum(stay_duration_before_merging) sum(stay_duration_after_merging)]];
        
        %% Cluster Time before-after AA-SPC on a cluster-level
        stay_duration_before_after_merging=[];
        for j=1:size(POP_sp_pair_ind_seq,1) % Sequence mereged clusters 
            dt_after_merging=ts_c(POP_sp_pair_ind_seq(j,2),2)-ts_c(POP_sp_pair_ind_seq(j,1),1);
            dt_before_merging=0;
            for k=POP_sp_pair_ind_seq(j,1):POP_sp_pair_ind_seq(j,2)
                dt_before_merging=dt_before_merging+(ts_c(k,2)-ts_c(k,1));
            end
            stay_duration_before_after_merging=[stay_duration_before_after_merging;[dt_before_merging dt_after_merging]];
        end
%         for j=1:length(single_POA_ind_after_seq_merged) % Unmerged clusters
%             dt=ts_c(single_POA_ind_after_seq_merged(j),2)-ts_c(single_POA_ind_after_seq_merged(j),1);
%             stay_duration_before_after_merging=[stay_duration_before_after_merging;[dt dt]];
%         end
        %% Discard Invalid Entries
        ind_test_temp=find(stay_duration_before_after_merging(:,2)-stay_duration_before_after_merging(:,1)<0);
        if(~isempty(ind_test_temp))
%             pid
%             stay_duration_before_after_merging(ind_test_temp,:)
            ind_v=find(stay_duration_before_after_merging(:,2)-stay_duration_before_after_merging(:,1)>=0);
            stay_duration_before_after_merging=stay_duration_before_after_merging(ind_v,:);
%             pause
        end
        POP_POA_cluster_level_stay_duration=[POP_POA_cluster_level_stay_duration;stay_duration_before_after_merging];
        %% Subject-level Avg. Cluster Time Diff before-after AA-SPC on a cluster-level
        stay_duration_before_after_merging=double(stay_duration_before_after_merging);
        percent_stay_duration_before_after_merging_diff=(stay_duration_before_after_merging(:,2)-stay_duration_before_after_merging(:,1))./stay_duration_before_after_merging(:,2)*100;
        POP_POA_subj_level_avg_stay_duration_diff=[POP_POA_subj_level_avg_stay_duration_diff;[pid mean(percent_stay_duration_before_after_merging_diff)]];
%         pause
    end
    %% List of subjects that has USEFUL (i.e., rectangle and radius find the same label) sp-pairs from the same bldg
    length(candidate_pids_sp_pair_D_C)
%     candidate_pids_sp_pair_D_C'
    
    %% Subject-level single-pair (i.e., NOT sequence) summary for POP, POA, POU --- each subject has at least 1 wlk session and each sp-pair has activity data
    size(POP_POA_POU_sp_pair_count_all_subj) % total 369 subjects 
    sum(POP_POA_POU_sp_pair_count_all_subj)
    
    %% Mean and SD of minute-level step counts during walk sessions across all subjects
    step_counts_wlk_all_subj=double(step_counts_wlk_all_subj);
    mean_wlk=mean(step_counts_wlk_all_subj); std_wlk=std(step_counts_wlk_all_subj); 
%     [mean_wlk std_wlk]
    range_wlk_generalized=[mean_wlk-std_wlk mean_wlk+std_wlk];
    
    %% Analysis of Cluster Count Decrease using AA-SPC
    POP_POA_sp_cluster_count_all_subj=POP_POA_POU_sp_cluster_count_all_subj; % new way of calculating CCD is based on POP, POA, and POU
    POP_POA_sp_cluster_count_all_subj=double(POP_POA_sp_cluster_count_all_subj);  % previous calculation of CCD didn't include POU
    size(POP_POA_sp_cluster_count_all_subj) % total 369 subjects
    sum(POP_POA_sp_cluster_count_all_subj)
    CCD=(POP_POA_sp_cluster_count_all_subj(:,2)-POP_POA_sp_cluster_count_all_subj(:,3))./POP_POA_sp_cluster_count_all_subj(:,2)*100;
    disp('Mean SD of CCD ....')
    [mean(CCD) std(CCD)]
    figure
    boxplot(CCD)
    grid on
    ylabel('% CCD')
%     xlabel('(d)')
    title('% Cluster Count Decrease (subject-level)')
    %% Write POP_POA_sp_cluster_count_all_subj into file
    outDir = [Dir 'outDir/AASPC_Evaluation_Summary_Unlabelled_NetHealth_v2/Personalized_walk_mean_std/'];
    if(~exist(outDir,'dir')); mkdir(outDir); end
    OutFile = [outDir 'AASPC_Cluster_Count_Before_After_Merging_CCD_per_subject' '.csv'];
    if(~isempty(POP_POA_sp_cluster_count_all_subj))
        fileID = fopen(OutFile,'w+');
        % pid,cluster_count(before),cluster_count(after),CCD (percent)
        fprintf(fileID,'pid,cluster_count(before),cluster_count(after),CCD (percent)\n');
        for j=1:size(POP_POA_sp_cluster_count_all_subj,1)
            fprintf(fileID,'%d,%f,%f,%f\n',POP_POA_sp_cluster_count_all_subj(j,:),CCD(j,:));
        end
        fclose(fileID);
    end
    
    %% Analysis of Cluster Time/Duration Increase using AA-SPC ---- cluster-level
    POP_POA_cluster_level_stay_duration=double(POP_POA_cluster_level_stay_duration);
    CTI=(POP_POA_cluster_level_stay_duration(:,2)-POP_POA_cluster_level_stay_duration(:,1))./POP_POA_cluster_level_stay_duration(:,2)*100;
    [h,p,ci,stats] = ttest(CTI/100)
    disp('Mean SD of CTI ...')
    [mean(CTI) std(CTI)]
    figure
    boxplot(CTI)
    grid on
    ylabel('% CTI')
%     xlabel('(d)')
    title('% Cluster Time Increase (cluster-level)')
    figure
    POP_POA_subj_level_avg_stay_duration_diff=double(POP_POA_subj_level_avg_stay_duration_diff);
    boxplot(POP_POA_subj_level_avg_stay_duration_diff(:,2))
    grid on
    ylabel('Avg. CTI (%)')
%     xlabel('(d)')
    title('Avg. Subject-level Cluster Time Increase (subject-level)')
    %% Write POP_POA_cluster_level_stay_duration into file
    OutFile = [outDir 'AASPC_Cluster_Time_Before_After_Merging_CTI_per_cluster' '.csv'];
    if(~isempty(POP_POA_cluster_level_stay_duration))
        fileID = fopen(OutFile,'w+');
        % cluster_time(before),cluster_time(after),CTI (percent)
        fprintf(fileID,'cluster_time(before),cluster_time(after),CTI (percent)\n');
        for j=1:size(POP_POA_cluster_level_stay_duration,1)
            fprintf(fileID,'%f,%f,%f\n',POP_POA_cluster_level_stay_duration(j,:),CTI(j,:));
        end
        fclose(fileID);
    end
    %% Write POP_POA_subj_level_avg_stay_duration_diff into file
    OutFile = [outDir 'AASPC_Avg_CTI_per_subject' '.csv'];
    if(~isempty(POP_POA_subj_level_avg_stay_duration_diff))
        fileID = fopen(OutFile,'w+');
        % pid,avg. CTI (percent)
        fprintf(fileID,'pid,avg. CTI (percent)\n');
        for j=1:size(POP_POA_subj_level_avg_stay_duration_diff,1)
            fprintf(fileID,'%d,%f\n',POP_POA_subj_level_avg_stay_duration_diff(j,:));
        end
        fclose(fileID);
    end
end
function [pid_wlk,ts_wlk,total_steps_wlk]=read_useful_wlk_sessions(File)
    pid_wlk=[]; ts_wlk=[]; 
    total_steps_wlk=[]; 
    if(~exist(File,'file')); disp('Can''t Read the walk sessions, labelled by Fitbit Smart Track Module xxx'); return
    else
        fid2 = fopen(File);
        % pid,st_ts(mili-sec),end_ts(mili-sec),total_steps
        B = textscan(fid2,'%f%f%f%f','delimiter',',','headerlines',1);
        fclose(fid2);
        [pid_wlk,st_ts_wlk,end_ts_wlk,total_steps_wlk]=B{1,:};
        ts_wlk = [st_ts_wlk end_ts_wlk];
    end
end