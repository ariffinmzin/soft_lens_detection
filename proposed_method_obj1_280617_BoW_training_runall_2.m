% this file is continuation from proposed_method_obj1_150617_BoW_runall.m
% this is for training the featureVector left and right using SVM
close all;
clear all;
clc;

global DIAGPATH
DIAGPATH = 'diagnostics';

load('soft.mat');
load('no.mat');

% read from excel file
[num,txt,raw] = xlsread('not_segmented_accurately.xlsx');
n=length(txt);
training_positive = zeros(n,500);
num_of_features = 1;

% count how many mat files in directory
files = dir('diagnostics/*.mat');
num_of_mat_files = numel(files);

% list all available mat files in a cell
list_of_filename_in_directory2 = cell(num_of_mat_files,3);

for k=1:num_of_mat_files
    
    parts = strsplit(files(k).name, '_');
    filename_in_directory = parts{3};
    parts2 = strsplit(filename_in_directory, '.');
    filename_in_directory2 = parts2{1};
    list_of_filename_in_directory2{k,1} = filename_in_directory2;
    
end % for k=1:num_of_mat_files

% load all mat files in not_segmented_accurately.xlsx - needed
for j=1:n

    try
        
        if (strcmp(txt{j,2},'left')==1) && (strcmp(txt{j,3},'')==1) % column 2 = 'left' and column 3 empty inside xlsx file 
    
            eval(sprintf('load(''featureVector_right_%s.mat'');',txt{j,1})); % true samples (purposely inversed)
            eval(sprintf('load(''featureVector_left_%s.mat'');',txt{j,1})); % negative samples
            
%             eval(sprintf('training_positive(num_of_features,:) = featureVector_right_%s(1,:);',txt{j,1}));
%             num_of_features = num_of_features + 1;
            
            for m=1:num_of_mat_files
                
                if (strcmp(list_of_filename_in_directory2{m,1},txt{j,1})==1)
                    
                    list_of_filename_in_directory2{m,2} = 1; % flag 1 for files that have been loaded from not_segmented_accurately.xlsx
                    list_of_filename_in_directory2{m,3} = 'left';
                    
                end
                
            end
            
        end
        
        if (strcmp(txt{j,3},'right')==1) && (strcmp(txt{j,2},'')==1)
            
            eval(sprintf('load(''featureVector_left_%s.mat'');',txt{j,1})); % true samples (purposely inversed)
            eval(sprintf('load(''featureVector_right_%s.mat'');',txt{j,1})); % negative samples
        
%             eval(sprintf('training_positive(num_of_features,:) = featureVector_left_%s(1,:);',txt{j,1}));
%             num_of_features = num_of_features + 1;
            
             for m=1:num_of_mat_files
                
                if (strcmp(list_of_filename_in_directory2{m,1},txt{j,1})==1)
                    
                    list_of_filename_in_directory2{m,2} = 1; % flag 1 for files that have been loaded from not_segmented_accurately.xlsx - needed
                    list_of_filename_in_directory2{m,3} = 'right';
                    
                end
                
             end

            
        end
        
%         num_of_features = num_of_features + 1;
%         eval(sprintf('training_positive(num_of_features,:) = featureVector_right_%s(1,:);',txt{j,1}));

        
    catch
            
        % purposely empty
 
    end
         
end % for j=1:n


%-----------------------------------------------------------------------------------------------------------------------


% load all mat files which segment accurately - SOFT LENS
for k=1:num_of_mat_files
    
    if (isempty(list_of_filename_in_directory2{k,2}) == 1) % empty at column 2
        
        eval(sprintf('load(''featureVector_left_%s.mat'');',list_of_filename_in_directory2{k,1})); % the subject name
        eval(sprintf('load(''featureVector_right_%s.mat'');',list_of_filename_in_directory2{k,1}));
        
    end
    
end

% for left (accurately segmented) - SOFT LENS
for s=1:num_of_mat_files
    
    if (isempty(list_of_filename_in_directory2{s,2}) == 1) && (isempty(list_of_filename_in_directory2{s,3}) == 1 && ismember(list_of_filename_in_directory2{s,1},soft) == 1)
        
        % 2135/1072 && ismember(list_of_filename_in_directory2{s,1},soft) == 1)% empty at column 2 and 3
        
        eval(sprintf('training_positive_left(%d,:) = featureVector_left_%s(1,:);',s,list_of_filename_in_directory2{s,1}));
%         eval(sprintf('training_positive_1 = featureVector_right_%s(1,:);',list_of_filename_in_directory2{k,1}));
              
    end
    
end

training_positive_left(~any(training_positive_left,2), :) = []; % remove any empty rows (zero elements)

% for s=1:num_of_mat_files

% for right (accurately segmented) - SOFT LENS
for s=1:num_of_mat_files
    
    if (isempty(list_of_filename_in_directory2{s,2}) == 1) && (isempty(list_of_filename_in_directory2{s,3}) == 1 && ismember(list_of_filename_in_directory2{s,1},soft) == 1)% empty at column 2 and 3
        
        eval(sprintf('training_positive_right(%d,:) = featureVector_right_%s(1,:);',s,list_of_filename_in_directory2{s,1}));
%         eval(sprintf('training_positive_1 = featureVector_right_%s(1,:);',list_of_filename_in_directory2{k,1}));
              
    end
    
end

training_positive_right(~any(training_positive_right,2), :) = []; % remove any empty rows (zero elements)

% for accurate right images from not_segmented_accurately.xlsx - SOFT LENS
for s=1:num_of_mat_files
    
    if (list_of_filename_in_directory2{s,2} == 1) & (strcmp(list_of_filename_in_directory2{s,3},'left') == 1 && ismember(list_of_filename_in_directory2{s,1},soft) == 1)
        
        eval(sprintf('training_positive_right_from_xlsx(%d,:) = featureVector_right_%s(1,:);',s,list_of_filename_in_directory2{s,1}));
%         eval(sprintf('training_positive_1 = featureVector_right_%s(1,:);',list_of_filename_in_directory2{k,1}));
              
    end
    
end

training_positive_right_from_xlsx(~any(training_positive_right_from_xlsx,2), :) = []; % remove any empty rows (zero elements)

% for accurate left images from not_segmented_accurately.xlsx - SOFT LENS
for s=1:num_of_mat_files
    
    if (list_of_filename_in_directory2{s,2} == 1) & (strcmp(list_of_filename_in_directory2{s,3},'right') == 1 && ismember(list_of_filename_in_directory2{s,1},soft) == 1)
        
        eval(sprintf('training_positive_left_from_xlsx(%d,:) = featureVector_left_%s(1,:);',s,list_of_filename_in_directory2{s,1}));
%         eval(sprintf('training_positive_1 = featureVector_right_%s(1,:);',list_of_filename_in_directory2{k,1}));
              
    end
    
end

training_positive_left_from_xlsx(~any(training_positive_left_from_xlsx,2), :) = []; % remove any empty rows (zero elements)

% concatenate all matrix (training_positive_left, training_positive_right,
% training_positive_right_from_xlsx, training_positive_left_from_xlsx)

training_data_positive = vertcat(training_positive_left, training_positive_right, training_positive_right_from_xlsx, training_positive_left_from_xlsx);



% for left images from not_segmented_accurately.xlsx (this is for negative
% samples - not perfectly segmented) - SOFT LENS
for s=1:num_of_mat_files
    
    if (list_of_filename_in_directory2{s,2} == 1) & (strcmp(list_of_filename_in_directory2{s,3},'left') == 1 && ismember(list_of_filename_in_directory2{s,1},soft) == 1)
        
        eval(sprintf('training_negative_left_from_xlsx(%d,:) = featureVector_left_%s(1,:);',s,list_of_filename_in_directory2{s,1}));
%         eval(sprintf('training_positive_1 = featureVector_right_%s(1,:);',list_of_filename_in_directory2{k,1}));
              
    end
    
end

training_negative_left_from_xlsx(~any(training_negative_left_from_xlsx,2), :) = []; % remove any empty rows (zero elements)

% for right images from not_segmented_accurately.xlsx (this is for negative
% samples - not perfectly segmented) - SOFT LENS
for s=1:num_of_mat_files
    
    if (list_of_filename_in_directory2{s,2} == 1) & (strcmp(list_of_filename_in_directory2{s,3},'right') == 1 && ismember(list_of_filename_in_directory2{s,1},soft) == 1)
        
        eval(sprintf('training_negative_right_from_xlsx(%d,:) = featureVector_right_%s(1,:);',s,list_of_filename_in_directory2{s,1}));
%         eval(sprintf('training_positive_1 = featureVector_right_%s(1,:);',list_of_filename_in_directory2{k,1}));
              
    end
    
end

training_negative_right_from_xlsx(~any(training_negative_right_from_xlsx,2), :) = []; % remove any empty rows (zero elements)

% need to have negative samples from the no lens images (8/12/17)
% need another column in list_of_filename_in_directory2 to indicate soft
% lens or without lens - this is covered by independent mat files
% (12/12/17)

% for left images from NO LENS - negative samples
for s=1:num_of_mat_files
    
    if (isempty(list_of_filename_in_directory2{s,2} == 1) && ismember(list_of_filename_in_directory2{s,1},no) == 1)
        
        eval(sprintf('training_negative_left(%d,:) = featureVector_left_%s(1,:);',s,list_of_filename_in_directory2{s,1}));
%         eval(sprintf('training_positive_1 = featureVector_right_%s(1,:);',list_of_filename_in_directory2{k,1}));
              
    end
    
end

training_negative_left(~any(training_negative_left,2), :) = []; % remove any empty rows (zero elements)

%2358 before remove any empty rows
%1064 after remove any empty rows

% for right images from NO LENS - negative samples
for s=1:num_of_mat_files
    
    if (isempty(list_of_filename_in_directory2{s,2} == 1) && ismember(list_of_filename_in_directory2{s,1},no) == 1)
        
        eval(sprintf('training_negative_right(%d,:) = featureVector_right_%s(1,:);',s,list_of_filename_in_directory2{s,1}));
%         eval(sprintf('training_positive_1 = featureVector_right_%s(1,:);',list_of_filename_in_directory2{k,1}));
              
    end
    
end

training_negative_right(~any(training_negative_right,2), :) = []; % remove any empty rows (zero elements)

%2358

training_data_negative = vertcat(training_negative_left_from_xlsx, training_negative_right_from_xlsx, training_negative_left, training_negative_right);

% negative 250 + 90

% 13/12/17 concatenate training data positive and training data negative
% X for data, Y for label

% 22/12/17
training_data_positive_negative = vertcat(training_data_positive, training_data_negative);
p_length = length(training_data_positive); % total number of y axis for positive data
n_length = length(training_data_negative); % total number of y axis for negative data
label = zeros(p_length+n_length,1);
label(1:p_length) = 1; % 1 for soft lens 0 for no lens

% training with SVM
%         SVMModel = fitcsvm(X,Y,'Standardize',true,'KernelFunction','RBF','KernelScale','auto');
        CVSVMModel = fitcsvm(training_data_positive_negative,label,'Standardize',true,'KernelFunction','RBF','KernelScale','auto','ClassNames',{'0','1'},'CrossVal','on');

%         [label,score] = predict(SVMModel,X)


% load('nd_data_training_LG4000.mat');
% n=length(nd_data_training_LG4000);
% i=1;
% soft_ndiris_training = cell(i,1);
% p=length(soft_ndiris_training); 
% k=1;
% 
% % load data
% 
% for j=1:n
%         
%     try
%     
%         if (strcmp(nd_data_training_LG4000{j,6},'Yes')==1)
%     
%             soft_ndiris_training{k,1} = nd_data_training_LG4000{j,1}; % filename
%     
%             eval(sprintf('load(''featureVector_left_%s.mat'');',soft_ndiris_training{k,1}));
%             
%             SVMModel = fitcsvm(X,Y,'Standardize',true,'KernelFunction','RBF','KernelScale','auto');
% 
%         end % if (strcmp(nd_data_training_LG4000{j,6},'Yes')==1)
%     catch
%     % empty
%     end
% 
% end % for j=1:n

