function [idx_clu,GMModel,clu_num,varargout] = sPCA_GMM_cluster_Color(data_in,varargin)

%varargin: bounds vector,K vector, number of time bins, pca vs spca vector,
%1 is spca

% define the control constant
if length(varargin) >= 9
    selector = varargin{9};
else
    selector = 0;
end

%if the user inputs arguments to run the sPCA also
if nargin > 2 && ~isempty(varargin{1})
    %extract the varargin contents
    bounds = varargin{1};
    K = varargin{2};
    t_bins = varargin{3};
    pca_vec = varargin{4};
    
    %% Generate data set for the sPCA and run it
    %define the number of groups
    g_num = size(bounds,2);
    %allocate memory to hold the aggregated data set 
    a_data = cell(g_num,1);
    
    %for each group of stimuli
    for stimg = 1:g_num
   
        %load the part of the trace corresponding to this stimulus    
        a_data{stimg} = data_in(:,bounds(1,stimg):bounds(2,stimg));
        
        %get the number of columns
        col_num = size(a_data{stimg},2);
        %for all the columns
        for col = 1:col_num
            
            %center the data by subtracting the mean of each column, and
            %normalizing each column by it's 2-norm
            a_data{stimg}(:,col) = (a_data{stimg}(:,col) - mean(a_data{stimg}(:,col)))...
                ./norm(a_data{stimg}(:,col)-mean(a_data{stimg}(:,col)));
        end
        %turn any NaNs and Infs into 0
        a_data{stimg}(isnan(a_data{stimg})) = 0;
        a_data{stimg}(isinf(a_data{stimg})) = 0;
        
    end
    %% Run the sPCA
    %Define parameters for the sPCA
    % K = [8 2 8 8 4 4 4];
    
    % K = [8 2 8 8];
    
    delta = 1000;
    % stop_var = -10;
    stop_var = -t_bins;
    maxSteps = 500;
    convergenceCriterion = 1e-9;
    verbose = false;
    
    %allocate memory for the PCs
    pcs = cell(g_num,5);
    
    %run the sPCA for each of the groups
    for stimg = 1:g_num
        %actual sPCA calculation
        [pcs{stimg,1},pcs{stimg,2},pcs{stimg,3},pcs{stimg,4},pcs{stimg,5}] = ...
            spca(a_data{stimg},[],K(stimg), delta, stop_var, maxSteps, convergenceCriterion, verbose);
    end
      
    %visualize spca results
    figure
    %for all the groups
    for stimg = 1:g_num
        subplot(ceil(sqrt(g_num)),ceil(sqrt(g_num)),stimg)
        
        if pca_vec(stimg) == 1
            imagesc(pcs{stimg,1})
        else
            imagesc(pcs{stimg,3}(:,1:K(stimg)))
        end
        title(strcat('PCs for'))
        ylabel('Time bins')
        xlabel('PC')
    end

    %load the pca results for export
    varargout{1} = pcs;
    %% Calculate the clustering vector using the sPCA results
    
    %get the length of the final feature vector
    ft_num = sum(K(1:g_num));
    f_bounds = cumsum([1,K]);
    %allocate memory for the composite vectors
    f_data = zeros(size(a_data{1},1),ft_num);
    
    %for all of the seeds
    for seeds = 1:size(f_data,1)
        %for all the groups of stimuli
        for stimg = 1:g_num
            if pca_vec(stimg) == 1
                %load the feature data matrix
                f_data(seeds,f_bounds(stimg):f_bounds(stimg+1)-1) = a_data{stimg}(seeds,:)*pcs{stimg,1};
            else
                %load the feature data matrix
                f_data(seeds,f_bounds(stimg):f_bounds(stimg+1)-1) = a_data{stimg}(seeds,:)*pcs{stimg,3}(:,1:K(stimg));
            end
        end
    end
    
    % if data was supplied for direct clustering combined with the spca
    % data
    if length(varargin) > 7 && ~isempty(varargin{8})
        % concatenate the matrices
        f_data = horzcat(f_data,varargin{8});
    end
    
    %standardize each column
    %for all the features
    for ft = 1:ft_num
        f_data(:,ft) = (f_data(:,ft) - mean(f_data(:,ft)))./std(f_data(:,ft));
    end
    figure
    % imagesc(f_data)
    % f_data(f_data>0.05) = 0.05;
    imagesc(log(abs(f_data)))
    title('Pre-clustering sparse PCs of the traces')
    xlabel('Time bins')
    ylabel('Individual traces')
    
    
    % if the selector is 2, exit
    if selector == 2
        idx_clu = [];
        GMModel = [];
        clu_num = [];
        varargout{3} = f_data;
        return
    end
else %if not with spca
    %just reassign the input data
    f_data = data_in;
    ft_num = size(f_data,2);
    %standardize each column
    %for all the features
    for ft = 1:ft_num
        f_data(:,ft) = (f_data(:,ft) - mean(f_data(:,ft)))./std(f_data(:,ft));
    end
end
%% Get the ideal cluster number based on the BIC


% %if the cluster number inputted is 0
% if clu_num == 0
%define the range of clu_num
clu_vec = varargin{6};

clu_ran = length(clu_vec);
%allocate memory to store the BIC for each clu_num
bic_vec = zeros(clu_ran,1);
%allocate memory for the models
model_vec = cell(clu_ran,1);
clu_count = 1;
%for all the clu num to evaluate
for clu = clu_vec
    fprintf(strjoin({'Current clu num:',num2str(clu_count),'of',num2str(length(clu_vec)),'\r\n'},'_'))
    % if the target cluster number is higher than the dimensionality, skip
    if clu >= size(f_data,1) || size(f_data,1) <= size(f_data,2)
        bic_vec(clu_count) = NaN;
        model_vec{clu_count} = [];
        clu_count = clu_count + 1;
        continue
    end
    %generate the model
    GM_temp = fitgmdist(f_data,clu,'CovarianceType','diagonal',...
        'RegularizationValue',1e-3,'Options',statset('MaxIter',500,...
        'Display','final'),'Replicates',varargin{7});
    
    bic_vec(clu_count) = GM_temp.BIC;
    model_vec{clu_count} = GM_temp;
    clu_count = clu_count + 1;
end

% Plot the BIC levels
figure
plot(clu_vec(2:end),2.*log10(abs(diff((bic_vec)))))
%send out the BIC levels
varargout{2} = bic_vec;
%assign the cluster number as the minimum of the BIC function
% check that there are values in the vector
if all(isnan(bic_vec))
    clu_num = [];
    GMModel = [];
    idx_clu = [];
else
    [~,clu_coord] = nanmin(bic_vec);
    clu_num = clu_vec(clu_coord);
    GMModel = model_vec{clu_coord};

    %cluster the data based on the model
    idx_clu = cluster(GMModel,f_data);
end