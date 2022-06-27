clc
clear all

path1 = '\traffic_data';
addpath(genpath([pwd,path1]))

filename = 'export4_light_complete_15min.csv';

raw = importdata(filename);

A = raw.textdata;
A(2:end,7) = num2cell(raw.data(:,1));