% Script for making Figure 7 and S9 (PLS regressions) in:

% Kardan, O., Stier, A. J., Cardenas-Iniguez, C., Schertz, K. E., Pruin, J. C., Deng, Y.,
% Chamberlain, T., Meredith, W. J.,  Zhang, X., Bowman, J. E., Lakhtakia, T., Tindel, L.,
% Avery, E. W., Yoo, K., Lin, Q., Chun, M. M., Berman, M. G., & Rosenberg, M. D.
% (accepted manuscript). Differences in the functional brain architecture of sustained attention
% and working memory in youth and adults. PLOS Biology 
% Data for this paper are at NDA (https://nda.nih.gov/) Study 1849 DOI: 10.15154/1528288

% Uses the PLS scripts from https://www.rotman-baycrest.on.ca/index.php?section=345
% Download plscmd and plsgui and place in Pls folder

% Uses bluewhitered from Nathan Childress (2022). 
% bluewhitered
% (https://www.mathworks.com/matlabcentral/fileexchange/4058-bluewhitered),
% MATLAB Central File Exchange.

% Uses violinplot from Holger Hoffmann (2022). 
% Violin Plot 
% (https://www.mathworks.com/matlabcentral/fileexchange/45134-violin-plot),
% MATLAB Central File Exchange. Retrieved December 5, 2022.

% questions? email Omid Kardan omidk@med.umich.edu

addpath(genpath('\Pls'));  % Pls folder
lowmotion =1; % set to 0 if running FD = .5 thresh for Figure S9
% Data are falttened Shen 268 correlation matrix from 0-back or 2-back
% stacked for each participant (5 columns of sub id and performance + 35778
% columns each represnting a pairwise connection. num Rwos equal sample size
% 
if lowmotion
    abcd0 = readtable('\shen_conn_matrices_abcd_stacked_subs_0bk.csv'); totn = 1545;
else
    abcd0 = readtable('\shen_conn_matrices_abcd_FD05_stacked_subs_0bk.csv');  totn = 3225;
end
hcp0 = readtable('\shen_conn_matrices_hcp_stacked_subs_0bk.csv');
if lowmotion
    abcd2 = readtable('\shen_conn_matrices_abcd_stacked_subs_2bk.csv'); totn = 1545;
else
    abcd2 = readtable('\shen_conn_matrices_abcd_FD05_stacked_subs_2bk.csv');  totn = 3225;
end
hcp2 = readtable('\shen_conn_matrices_hcp_stacked_subs_2bk.csv');
%%
%%%%%%%%% SA
[a1,b] = find(isnan(abcd0{:,6:35783})); %to remove participants with NaN connections because PLs won't run
[a2,b] = find(isnan(hcp0{:,6:35783})); %
goodabcdsubs = setdiff([1:totn],a1); % should be 1502 for low-motion and 3135 for high-motion PLS
goodhcpsubs = setdiff([1:754],a2); % should stay 754

rng('shuffle');
for k=1:200
    k_subs = randperm(length(goodabcdsubs),length(goodhcpsubs));
    goodabcdsubs_k = goodabcdsubs(k_subs);
    groups{1} = [goodabcdsubs_k];   % children
    groups{2} = [goodhcpsubs];  % adults
    
    nRuns = 1;Runs=[1];
    nParcels = 268; 
    goodps = 1:268;
    inds = find(tril(ones(length(goodps)),-1)==1);
    
    datamat_lst = cell(length(groups),1); delaprimes =[abcd0.Acc(goodabcdsubs_k); hcp0.Acc(goodhcpsubs)]; % accuracy in x-back
    
    for g = 1:numel(groups)
        for c = 1:nRuns
            
            if g==1
                datamat_lst{g} = abcd0{goodabcdsubs_k,6:35783};
            end
            if g==2
                datamat_lst{g} = hcp0{goodhcpsubs,6:35783};
            end
            g
            
        end
    end
    
    
    num_subj = [length(groups{1})  length(groups{2})];
    num_cond = nRuns;
    option.method = 3; % behavioral PLS3
    option.num_boot = 500;
    option.num_perm = 500;
    option.meancentering_type=[2];
    option.cormode = 0; %	0. Pearson correlation 2. covaraince 4. cosine angle 6. dot product
    option.stacked_behavdata = [delaprimes]
    result = pls_analysis(datamat_lst, num_subj, num_cond, option);
    save(['result_',num2str(k),'.mat'],'result');

%%%%%%%%% WM
[a1,b] = find(isnan(abcd2{:,6:35783})); %to remove participants with NaN connections because PLs won't run
[a2,b] = find(isnan(hcp2{:,6:35783})); %
goodabcdsubs = setdiff([1:totn],a1); % should be 1502 for low-motion and 3135 for high-motion PLS
goodhcpsubs = setdiff([1:length(goodhcpsubs)],a2); % should stay 754

    goodabcdsubs_k = goodabcdsubs(k_subs);
    groups{1} = [goodabcdsubs_k];   % children
    groups{2} = [goodhcpsubs];  % adults
    
    nRuns = 1;Runs=[1];
    nParcels = 268;
    goodps = 1:268;
    inds = find(tril(ones(length(goodps)),-1)==1);
    
    datamat_lst = cell(length(groups),1); delaprimes =[abcd2.Acc(goodabcdsubs_k); hcp2.Acc(goodhcpsubs)]; % accuracy in 2-back
    
    for g = 1:numel(groups)
        for c = 1:nRuns
            
            if g==1
                datamat_lst{g} = abcd2{goodabcdsubs_k,6:35783};
            end
            if g==2
                datamat_lst{g} = hcp2{goodhcpsubs,6:35783};
            end
            g
            
        end
    end
    
    
    num_subj = [length(groups{1})  length(groups{2})];
    num_cond = nRuns;
    option.method = 3; %behavioral PLS3
    option.num_boot = 500;
    option.num_perm = 500;
    option.meancentering_type=[2];
    option.cormode = 0; %	0. Pearson correlation 2. covaraince 4. cosine angle 6. dot product
    option.stacked_behavdata = [delaprimes]
    result = pls_analysis(datamat_lst, num_subj, num_cond, option);
    save(['resultwm_',num2str(k),'.mat'],'result');
end

%% SA aggregate plot
clear all
    goodps = 1:268;
    inds = find(tril(ones(length(goodps)),-1)==1);
addpath(genpath('\bluewhitered'));
addpath(genpath('\Violinplot-Matlab-master'));
SAsprob =[];
SAdistrib_11 = [];
SAdistrib_21 = [];
SAdistrib_12 = [];
SAdistrib_22 = [];
SAs = [];
SAcompare_u_1 = [];
SAcompare_u_2 = [];

for i=1:200
    load(['result_',num2str(i),'.mat'])
    SAsprob = [SAsprob result.perm_result.sprob];
    
    SAdistrib_11 = [SAdistrib_11 squeeze(result.boot_result.distrib(1,1,:))];
    SAdistrib_21 = [SAdistrib_21 squeeze(result.boot_result.distrib(2,1,:))];
    
    SAdistrib_12 = [SAdistrib_12 squeeze(result.boot_result.distrib(1,2,:))];
    SAdistrib_22 = [SAdistrib_22 squeeze(result.boot_result.distrib(2,2,:))];
    
    SAs = [SAs result.s];
    SAcompare_u_1 = [SAcompare_u_1 result.boot_result.compare_u(:,1)];
    
    SAcompare_u_2 = [SAcompare_u_2 result.boot_result.compare_u(:,2)];
     
end

    ps = SAsprob
    figure;
    subplot(2,2,1)
    % first LV
    lv=1;
    flipped = find(mean(SAdistrib_11<0) & mean(SAdistrib_21<0)); % check for 180 degrees rotation
    SAdistrib_11(:,flipped) = (-1)*SAdistrib_11(:,flipped);
    SAdistrib_21(:,flipped) = (-1)*SAdistrib_21(:,flipped);
    SAcompare_u_1(:,flipped) = (-1)*SAcompare_u_1(:,flipped);
    
    violinplot([reshape(SAdistrib_11,numel(SAdistrib_11),1) reshape(SAdistrib_21,numel(SAdistrib_21),1)]);
    xlim([.5,2.5]);ylim([-.65,1]);
    
    cbcovs = mean(SAs,2).^2./sum(mean(SAs,2).^2);
    title(['LV',num2str(lv),'  \sigma_{XY} = ',num2str(cbcovs(lv)),'  p = ',num2str(ps(lv))]);
    ylabel('Correlation with 0-back performance','FontSize',18);
    set(gca,'XTick',1:2,...
        'Xticklabel',{' ',' '},'XAxisLocation','origin',...
        'XtickLabelRotation',45,'FontSize',14);
    hold off
    axis square
    subplot(2,2,2)
    SAcompare_u_1 = mean(SAcompare_u_1,2);
    brainZs1 = SAcompare_u_1(abs(SAcompare_u_1)>3);
    temp1=zeros(268,268);temp1(inds(abs(SAcompare_u_1)>3))=brainZs1;
    max(max(temp1))
    min(min(temp1))
    
    imagesc(temp1,[-8.8,7.73]); axis square
    
    colormap(bluewhitered), colorbar;
    xnames ={'R-prefrontal ' 'R-motor ' '\rightarrow' 'R-parietal ' 'R-temporal ' 'R-occipital ' 'R-limbic '...
        'R-cerebellum ' 'R-subcortex ' '-\rightarrow' 'L-prefrontal ' 'L-motor ' '\rightarrow' 'L-parietal '...
        'L-temporal ' 'L-occipital ' 'L-limbic ' 'L-cerebellum ' 'L-subcortex '  '-\rightarrow'};

    netsizes1 = [0 22 11 4 13 21 11 17 20 9 5 24 10 3 14 18 14 19 21 8];
    netsizes = [22 11 4 13 21 11 17 20 9 5 24 10 3 14 18 14 19 21 8 4];
    set(gca,'XTick',round(cumsum(netsizes1)+ netsizes/2),'Xticklabel',xnames,...
        'XtickLabelRotation',90,'FontSize',13);
    set(gca,'YTick',round(cumsum(netsizes1)+ netsizes/2),'Yticklabel',xnames,...
        'YtickLabelRotation',0,'FontSize',13,'TickLength',[0,0]);
    
    makans = cumsum(netsizes);
    for j=1:length(netsizes)
        line([makans(j),makans(j)],[makans(j) ,333],'Color','black');
    end
    for j=1:length(netsizes)
        line([0 ,makans(j)],[makans(j),makans(j)],'Color','black');
    end
    
    % second LV
    subplot(2,2,3)
    lv=2;
    
    flipped = find(mean(SAdistrib_12<0) & mean(SAdistrib_22>0)); % check for 180 degrees rotation
    SAdistrib_12(:,flipped) = (-1)*SAdistrib_12(:,flipped);
    SAdistrib_22(:,flipped) = (-1)*SAdistrib_22(:,flipped);
    SAcompare_u_2(:,flipped) = (-1)*SAcompare_u_2(:,flipped);
    
    violinplot([reshape(SAdistrib_12,numel(SAdistrib_12),1) reshape(SAdistrib_22,numel(SAdistrib_22),1)]);
    xlim([.5,2.5]);ylim([-1,1]);
    
    cbcovs = mean(SAs,2).^2./sum(mean(SAs,2).^2);
    title(['LV',num2str(lv),'  \sigma_{XY} = ',num2str(cbcovs(lv)),'  p = ',num2str(ps(lv))]);
    ylabel('Correlation with 0-back performance','FontSize',18);
    set(gca,'XTick',1:2,...
        'Xticklabel',{' ',' '},'XAxisLocation','origin',...
        'XtickLabelRotation',45,'FontSize',14);
    axis square
    hold off
    subplot(2,2,4)
    
        SAcompare_u_2 = mean(SAcompare_u_2,2);
    brainZs2 = SAcompare_u_2(abs(SAcompare_u_2)>3);
    temp2=zeros(268,268);temp2(inds(abs(SAcompare_u_2)>3))=brainZs2;
    max(max(temp2))
    min(min(temp2))
    
    imagesc(temp2,[-8.8,7.73]);
    
    axis square
    colormap(bluewhitered), colorbar;

    set(gca,'XTick',round(cumsum(netsizes1)+ netsizes/2),...
    'Xticklabel',xnames,'XtickLabelRotation',90,'FontSize',13);
    set(gca,'YTick',round(cumsum(netsizes1)+ netsizes/2),'Yticklabel',xnames,...
        'YtickLabelRotation',0,'FontSize',13,'TickLength',[0,0]);
    
    makans = cumsum(netsizes);
    for j=1:length(netsizes)
        line([makans(j),makans(j)],[makans(j) ,333],'Color','black');
    end
    for j=1:length(netsizes)
        line([0 ,makans(j)],[makans(j),makans(j)],'Color','black');
    end

    %%
    
%% WM aggregate plots
clear all
    goodps = 1:268;
    inds = find(tril(ones(length(goodps)),-1)==1);
addpath(genpath('\bluewhitered'));
addpath(genpath('\Violinplot-Matlab-master'));

WMsprob = [];
WMdistrib_11 = [];
WMdistrib_21 = [];
WMdistrib_12 = [];
WMdistrib_22 = [];
WMs = [];
WMcompare_u_1 = [];
WMcompare_u_2 = [];

for i=1:200
        
    load(['resultwm_',num2str(i),'.mat'])
    WMsprob = [WMsprob result.perm_result.sprob];
    
    WMdistrib_11 = [WMdistrib_11 squeeze(result.boot_result.distrib(1,1,:))];
    WMdistrib_21 = [WMdistrib_21 squeeze(result.boot_result.distrib(2,1,:))];
    
    WMdistrib_12 = [WMdistrib_12 squeeze(result.boot_result.distrib(1,2,:))];
    WMdistrib_22 = [WMdistrib_22 squeeze(result.boot_result.distrib(2,2,:))];
    
    WMs = [WMs result.s];
    WMcompare_u_1 = [WMcompare_u_1 result.boot_result.compare_u(:,1)];
    
    WMcompare_u_2 = [WMcompare_u_2 result.boot_result.compare_u(:,2)];
    
end

    ps = WMsprob
    figure;
    subplot(2,2,1)
    % first LV
    lv=1;
    flipped = find(mean(WMdistrib_11<0) & mean(WMdistrib_21<0)); % check for 180 degrees rotation

    WMdistrib_11(:,flipped) = (-1)*WMdistrib_11(:,flipped);
    WMdistrib_21(:,flipped) = (-1)*WMdistrib_21(:,flipped);
    WMcompare_u_1(:,flipped) = (-1)*WMcompare_u_1(:,flipped);
    
    violinplot([reshape(WMdistrib_11,numel(WMdistrib_11),1) reshape(WMdistrib_21,numel(WMdistrib_21),1)]);
    xlim([.5,2.5]);ylim([-.65,1]);
    
    cbcovs = mean(WMs,2).^2./sum(mean(WMs,2).^2);
    title(['LV',num2str(lv),'  \sigma_{XY} = ',num2str(cbcovs(lv)),'  p = ',num2str(ps(lv))]);
    ylabel('Correlation with 2-back performance','FontSize',18);
    set(gca,'XTick',1:2,...
        'Xticklabel',{' ',' '},'XAxisLocation','origin',...
        'XtickLabelRotation',45,'FontSize',14);
    hold off
    axis square
    subplot(2,2,2)
    WMcompare_u_1 = mean(WMcompare_u_1,2);
    brainZs1 = WMcompare_u_1(abs(WMcompare_u_1)>3);
    temp1=zeros(268,268);temp1(inds(abs(WMcompare_u_1)>3))=brainZs1;
    max(max(temp1))
    min(min(temp1))
    imagesc(temp1,[-9,7.64]); axis square
    colormap(bluewhitered), colorbar;
    xnames ={'R-prefrontal ' 'R-motor ' '\rightarrow' 'R-parietal ' 'R-temporal ' 'R-occipital ' 'R-limbic '...
        'R-cerebellum ' 'R-subcortex ' '-\rightarrow' 'L-prefrontal ' 'L-motor ' '\rightarrow' 'L-parietal '...
        'L-temporal ' 'L-occipital ' 'L-limbic ' 'L-cerebellum ' 'L-subcortex '  '-\rightarrow'};

    netsizes1 = [0 22 11 4 13 21 11 17 20 9 5 24 10 3 14 18 14 19 21 8];
    netsizes = [22 11 4 13 21 11 17 20 9 5 24 10 3 14 18 14 19 21 8 4];
    set(gca,'XTick',round(cumsum(netsizes1)+ netsizes/2),'Xticklabel',xnames,...
        'XtickLabelRotation',90,'FontSize',13);
    set(gca,'YTick',round(cumsum(netsizes1)+ netsizes/2),'Yticklabel',xnames,...
        'YtickLabelRotation',0,'FontSize',13);
    
    makans = cumsum(netsizes);
    for j=1:length(netsizes)
        line([makans(j),makans(j)],[makans(j) ,333],'Color','black');
    end
    for j=1:length(netsizes)
        line([0 ,makans(j)],[makans(j),makans(j)],'Color','black');
    end
    
    % second LV
    subplot(2,2,3)
    lv=2;
    flipped = find(mean(WMdistrib_12<0) & mean(WMdistrib_22>0)); % check for 180 degrees rotation
    WMdistrib_12(:,flipped) = (-1)*WMdistrib_12(:,flipped);
    WMdistrib_22(:,flipped) = (-1)*WMdistrib_22(:,flipped);
    WMcompare_u_2(:,flipped) = (-1)*WMcompare_u_2(:,flipped);
    
    violinplot([reshape(WMdistrib_12,numel(WMdistrib_12),1) reshape(WMdistrib_22,numel(WMdistrib_22),1)]);
    xlim([.5,2.5]);ylim([-1,1]);
    
    cbcovs = mean(WMs,2).^2./sum(mean(WMs,2).^2);
    title(['LV',num2str(lv),'  \sigma_{XY} = ',num2str(cbcovs(lv)),'  p = ',num2str(ps(lv))]);
    ylabel('Correlation with 2-back performance','FontSize',18);
    set(gca,'XTick',1:2,...
        'Xticklabel',{' ',' '},'XAxisLocation','origin',...
        'XtickLabelRotation',45,'FontSize',14);
    axis square
    hold off
    subplot(2,2,4)
    
        WMcompare_u_2 = mean(WMcompare_u_2,2);
    brainZs2 = WMcompare_u_2(abs(WMcompare_u_2)>3);
    temp2=zeros(268,268);temp2(inds(abs(WMcompare_u_2)>3))=brainZs2;
    max(max(temp2))
    min(min(temp2))
    
    imagesc(temp2,[-8.8,7.73]);
    
    axis square
    colormap(bluewhitered), colorbar;

    set(gca,'XTick',round(cumsum(netsizes1)+ netsizes/2),'Xticklabel',xnames,...
        'XtickLabelRotation',90,'FontSize',13);
    set(gca,'YTick',round(cumsum(netsizes1)+ netsizes/2),'Yticklabel',xnames,...
        'YtickLabelRotation',0,'FontSize',13);
    
    makans = cumsum(netsizes);
    for j=1:length(netsizes)
        line([makans(j),makans(j)],[makans(j) ,333],'Color','black');
    end
    for j=1:length(netsizes)
        line([0 ,makans(j)],[makans(j),makans(j)],'Color','black');
    end

