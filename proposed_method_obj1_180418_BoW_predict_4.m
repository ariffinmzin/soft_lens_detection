close all;
clear all;
clc;

global DIAGPATH
DIAGPATH = 'diagnostics';

load('SVMModel_221217.mat');
CompactSVMModel = CVSVMModel.Trained{1}; % Extract trained, compact classifier

% load training data

%     load('nd_data_training_LG4000.mat');
%     n=length(nd_data_training_LG4000);
%     i=1;
%     soft_ndiris_training = cell(i,7);
%     p=length(soft_ndiris_training); 
%     k=1;
%     
%     right_training = 0;
%     left_training = 0;
%     
%     for j=1:n
%       
%         % try
%            
%          if (strcmp(nd_data_training_LG4000{j,6},'Yes')==1)
%             
%             soft_ndiris_training{k,1} = nd_data_training_LG4000{j,1}; % filename
%             soft_ndiris_training{k,2} = nd_data_training_LG4000{j,7}; % pupilx
%             soft_ndiris_training{k,3} = nd_data_training_LG4000{j,8}; % pupily
%             soft_ndiris_training{k,4} = nd_data_training_LG4000{j,9}; % pupilr
%             soft_ndiris_training{k,5} = nd_data_training_LG4000{j,10}; % irisx
%             soft_ndiris_training{k,6} = nd_data_training_LG4000{j,11}; % irisy
%             soft_ndiris_training{k,7} = nd_data_training_LG4000{j,12}; % irisr
%             
%             display(soft_ndiris_training{k,1});
%             
%             try
%                 
%                 eval(sprintf('load(''featureVector_right_%s.mat'');',soft_ndiris_training{k,1})); 
%                 right_training = right_training + 1; % to count how many successful feature extracted
%                 % predict
%                 % [label,score] = predict(CVSVMModel,strcat(featureVector_right_,soft_ndiris_training{k,1}));
%             catch
%                 
%             end
%             
%             try
%                 
%                 eval(sprintf('load(''featureVector_left_%s.mat'');',soft_ndiris_training{k,1}));
%                 left_training = left_training + 1;
%                 
%             catch
%                 
%             end
%             
%          end % if
%          
%     end
         
%% 
  
% load testing data - with soft lens
    load('nd_data_testing_LG4000.mat');
    n=length(nd_data_testing_LG4000);
    i=1;
    soft_ndiris_testing = cell(i,7);
    p=length(soft_ndiris_testing); 
    k=1;
    
    right_testing_soft = 0;
    left_testing_soft = 0;
%     result_right = cell(300,2);
    
    for j=1:n
      
        % try
           
         if (strcmp(nd_data_testing_LG4000{j,6},'Yes')==1)
            
            soft_ndiris_testing{k,1} = nd_data_testing_LG4000{j,1}; % filename
            soft_ndiris_testing{k,2} = nd_data_testing_LG4000{j,7}; % pupilx
            soft_ndiris_testing{k,3} = nd_data_testing_LG4000{j,8}; % pupily
            soft_ndiris_testing{k,4} = nd_data_testing_LG4000{j,9}; % pupilr
            soft_ndiris_testing{k,5} = nd_data_testing_LG4000{j,10}; % irisx
            soft_ndiris_testing{k,6} = nd_data_testing_LG4000{j,11}; % irisy
            soft_ndiris_testing{k,7} = nd_data_testing_LG4000{j,12}; % irisr
            
            display(soft_ndiris_testing{k,1});
            
            % right sclera
            try
                
                eval(sprintf('load(''featureVector_right_%s.mat'');',soft_ndiris_testing{k,1})); 
                right_testing_soft = right_testing_soft + 1; % to count how many successful feature extracted
                % combine all features in a vector 
                eval(sprintf('testing_positive_right(%d,:) = featureVector_right_%s(1,:);',j,soft_ndiris_testing{k,1}));
                

                
                
            catch
                
            end
            
            
            % soft_ndiris_testing = cell(right_testing,7);
            
            % left sclera
            try
                
                eval(sprintf('load(''featureVector_left_%s.mat'');',soft_ndiris_testing{k,1}));
                left_testing_soft = left_testing_soft + 1;
                % combine all features in a vector 
                eval(sprintf('testing_positive_left(%d,:) = featureVector_left_%s(1,:);',j,soft_ndiris_testing{k,1}));
                
            catch
                
            end
         end % if
         
         
        % catch
            
        % end
         
    end
    
    % predict testing data right - with soft lens
    
    testing_positive_right(~any(testing_positive_right,2), :) = []; % remove any empty rows (zero elements)
    
    eval(sprintf('[label_testing_positive_right, score_testing_positive_right] = predict(CompactSVMModel,testing_positive_right);'));
    
    % predict testing data left - with soft lens 
    
    testing_positive_left(~any(testing_positive_left,2), :) = []; % remove any empty rows (zero elements)
    
    eval(sprintf('[label_testing_positive_left, score_testing_positive_left] = predict(CompactSVMModel,testing_positive_left);'));
    
    % count how many TP
    count_TP_right = sum(label_testing_positive_right(:) == 1);
    
    % count how many TP
    count_TP_left = sum(label_testing_positive_left(:) == 1);
    
    %% 
    
% load testing data - without soft lens
    load('nd_data_testing_LG4000.mat');
    n=length(nd_data_testing_LG4000);
    i=1;
    no_ndiris_testing = cell(i,7);
    p=length(no_ndiris_testing); 
    k=1;
    
    right_testing_no = 0;
    left_testing_no = 0;
%     result_right = cell(300,2);
    
    for j=1:n
      
        % try
           
         if (strcmp(nd_data_testing_LG4000{j,6},'No')==1)
            
            no_ndiris_testing{k,1} = nd_data_testing_LG4000{j,1}; % filename
            no_ndiris_testing{k,2} = nd_data_testing_LG4000{j,7}; % pupilx
            no_ndiris_testing{k,3} = nd_data_testing_LG4000{j,8}; % pupily
            no_ndiris_testing{k,4} = nd_data_testing_LG4000{j,9}; % pupilr
            no_ndiris_testing{k,5} = nd_data_testing_LG4000{j,10}; % irisx
            no_ndiris_testing{k,6} = nd_data_testing_LG4000{j,11}; % irisy
            no_ndiris_testing{k,7} = nd_data_testing_LG4000{j,12}; % irisr
            
            display(no_ndiris_testing{k,1});
            
            % right sclera
            try
                
                eval(sprintf('load(''featureVector_right_%s.mat'');',no_ndiris_testing{k,1})); 
                right_testing_no = right_testing_no + 1; % to count how many successful feature extracted
                % combine all features in a vector 
                eval(sprintf('testing_negative_right(%d,:) = featureVector_right_%s(1,:);',j,no_ndiris_testing{k,1}));     
                
            catch
                
            end
            
            
            % soft_ndiris_testing = cell(right_testing,7);
            
            % left sclera
            try
                
                eval(sprintf('load(''featureVector_left_%s.mat'');',no_ndiris_testing{k,1}));
                left_testing_no = left_testing_no + 1;
                % combine all features in a vector 
                eval(sprintf('testing_negative_left(%d,:) = featureVector_left_%s(1,:);',j,no_ndiris_testing{k,1}));
                
            catch
                
            end
         end % if
         
         
        % catch
            
        % end
         
    end
    
    % predict testing data right - without soft lens
    
    testing_negative_right(~any(testing_negative_right,2), :) = []; % remove any empty rows (zero elements)
    
    eval(sprintf('[label_testing_negative_right, score_testing_negative_right] = predict(CompactSVMModel,testing_negative_right);'));
    
    % predict testing data left - without soft lens 
    
    testing_negative_left(~any(testing_negative_left,2), :) = []; % remove any empty rows (zero elements)
    
    eval(sprintf('[label_testing_negative_left, score_testing_negative_left] = predict(CompactSVMModel,testing_negative_left);'));
    
    % count how many TN
    count_TN_right = sum(label_testing_negative_right(:) == 1);
    
    % count how many TN
    count_TN_left = sum(label_testing_negative_left(:) == 1);
    
    

