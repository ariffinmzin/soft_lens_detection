% created on 21/05/17 - to train the line intensity profie (segmented lens
% boundary pixels) with SVM

close all;
clear all;
clc;

global DIAGPATH
DIAGPATH = 'diagnostics';

% load('nd_data_training_LG4000.mat');
% n=length(nd_data_training_LG4000);
% i=1;
% soft_ndiris_training = cell(i,7);
% p=length(soft_ndiris_training); 
% k=1;
% 
% load(strcat('left_lens_boundary_xy_','.mat'));

% load soft data from directory soft
files = dir('diagnostics/soft/*.mat');

num_files = numel(files);
xname=zeros(num_files,1);
for w=1:num_files
    
    xname=getfield(files,{w,1},'name');
%     load(strcat(xname,'.mat'));
    load(xname);

end

% load no data from directory no
% files = dir('diagnostics/no/*.mat');
% 
% num_files = numel(files);
% xname=zeros(num_files,1);
% for w=1:num_files
%     
%     xname=getfield(files,{w,1},'name');
%     load(strcat(xname,'.mat'));
% 
% end


