%analysis

clc
clear
close all

% Add helper functions
addpath(genpath(fullfile(pwd, 'helperFunctions')));

% Load the data
load data

% Boundaries for analyzing block halves or quantiles
bounds2 = [0, 25, 50];
bounds5 = [0, 10, 20, 30, 40, 50];

%% Loop over all subjects and compute basic quantities
for sub=1:length(data)
    
    % Simplify data
    stim = data{sub}.stim;
    resp = data{sub}.resp;
    conf = data{sub}.conf;
    correct = data{sub}.correct;
    condition = data{sub}.condition;
    rt_raw = data{sub}.rt;
    angle = data{sub}.angle;
    trial_in_block = data{sub}.trial_in_block;
        
    % Compute basic performance metrics
    for cond=1:2
        [dprime(sub,cond),c(sub,cond,:), log_beta(sub,cond,:)] = computeSDTcriteria(...
            stim(condition==cond), resp(condition==cond), conf(condition==cond), 4);
        confidence(sub,cond) = mean(conf(condition==cond));
        accuracy(sub,cond) = mean(correct(condition==cond));
        rt(sub,cond) = mean(rt_raw(condition==cond));
        
        % Split-half analyses
        for part=1:length(bounds2)-1
            conf_halves(cond,part,sub) = mean(conf(condition==cond & trial_in_block>bounds2(part) &...
                trial_in_block<=bounds2(part+1)));
        end
        
        % Block quantiles
        for part=1:length(bounds5)-1
            conf_quantiles(cond,part,sub) = mean(conf(condition==cond & trial_in_block>bounds5(part) &...
                trial_in_block<=bounds5(part+1)));
        end
    end
    phi(sub) = corr(conf', correct');
end
% Remove Inf and -Inf from c
c(c==-Inf) = NaN; c(c==Inf) = NaN;


%% Compare the criteria between the two conditions
fprintf('\n--------- Compare the criteria between the two conditions ---------');
t(:,1,:) = c(:,1,:)*6; %threshold in degrees
t(:,2,:) = c(:,2,:)*18;
average_t0 = mean(t(:,:,4))
perform_ttest(t(:,:,4), 'Compare t_0 to 0');
perform_ttest(t(:,2,4)-t(:,1,4), 't_0 comparison');
r = squeeze(t(:,2,:)./t(:,1,:));
for crit_loc=[1:3,5:7]
    [p_nonparametric(crit_loc)] = signrank(r(:,crit_loc)-3);
end
p_nonparametric
perform_ttest(squeeze(t(:,2,:))-3*squeeze(t(:,1,:)), 'criterion');


%% Basic comparisons between the two conditions
fprintf('\n--------- Basic comparisons between the two conditions ---------');
mean_dprime = mean(dprime)
perform_ttest(dprime(:,2)-dprime(:,1), 'd''');
mean_confidence = mean(confidence)
perform_ttest(confidence(:,2)-confidence(:,1), 'confidence');
confEffect_pos_neg = [sum(diff(confidence,1,2)'>0), sum(diff(confidence,1,2)'<0)]
mean_rt = squeeze(mean(rt))
perform_ttest(squeeze(rt(:,2,:))-squeeze(rt(:,1,:)), 'RT');
perform_corr(diff(confidence,1,2), diff(dprime,1,2), 'Relationship b/n confidence effect and d'' effect', 1);
perform_corr(diff(confidence,1,2), diff(rt,1,2), 'Relationship b/n confidence effect and RT effect', 1);


%% Model fitting results
fprintf('\n--------- Model results ---------');
load modeling/model_fits
expansion_coef = params(7,:);
attraction_coef = sqrt(3./(expansion_coef.^2));
mean_min_max_SD_attraction = [mean(attraction_coef), min(attraction_coef), max(attraction_coef), std(attraction_coef)]
percent_movement = [mean_min_max_SD_attraction(1)-1, sqrt(3)-1]*100
perform_ttest(attraction_coef-1, 'Attraction');
perform_corr(attraction_coef', diff(confidence,1,2), 'Relationship b/n attraction and confidence effect', 1);


%% Evolution of the effects over the course of the block
fprintf('\n--------- Separating the block into different subparts ---------');
mean_conf_halves = mean(conf_halves,3)
conf_halves_diff = squeeze(diff(conf_halves))';
conf_diff_firstHalf_secondHalf = mean(conf_halves_diff)
perform_ttest(conf_halves_diff(:,1)-conf_halves_diff(:,2), 'Second-first half conf for Low Var cond');
mean_conf_quantiles = diff(mean(conf_quantiles,3))


%% Plot Figures 3-5 from paper
% % Figure 3
% plot_criteria(t);
% 
% % Figure 4
% figure('Color','w', 'DefaultAxesFontSize',16);
% plot_2bars(dprime, 'd''', [.3,1.41], 1);
% plot_2bars(confidence, 'Confidence', [1.5,4], 2);
% plot_2bars(rt, 'RT (sec)', [.5,1.8],3);
% plot_scatterplot(diff(confidence,1,2), diff(dprime,1,2), 4, 'conf_{HighVar}-conf_{LowVar}', 'd''_{HighVar}-d''_{LowVar}');
% plot_scatterplot(diff(confidence,1,2), diff(rt,1,2), 6, 'conf_{HighVar}-conf_{LowVar}', 'RT_{HighVar}-RT_{LowVar}');
% 
% % Figure 5
% plot_1bar(attraction_coef, 'criterion attraction (\alpha)', [.8,1.7]);
% plot_scatterplot(diff(confidence,1,2), attraction_coef', NaN, 'conf_{HighVar}-conf_{LowVar}', 'criterion attraction (\alpha)')