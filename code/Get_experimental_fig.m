% compare different result
% @author Wang Yangyang
% 

addpath(genpath('tools'));
dbstop if error;
clear;
clc;
disp('show result');

names = {
    { 'Lena';     38000; };
    { 'Baboon';   13000; };
    { 'Airplane'; 51000; }; 
    { 'Boat';     26000; };
    { 'Peppers';  31000; };
    { 'Lake';     26000; };
    { 'Barbara';  31000; };
    { 'Elaine';   24000; };
};

mtds = {
%     {'Proposed'; '2018'; ':o';  'k'};
    {'Proposed'; '2019'; '-o';  'k'};
    {'PVO';  '2013'; ':o'; 'r'}; 
    {'IPVO'; '2013'; ':o';  'b'};
    {'PVOK'; '2014'; ':o';  'g'};
%     {'PairwisePVO'; '2018'; '-o';  'b'};
%     {'PPVO'; '2016'; ':o';  'b'};
%     {'APPVO'; '2018'; '-o';  'k'};
%     {'SachnevEC'; '2009'; '-o';  'c'};
%     {'M2WPVO'; '2018'; '-v';  'y'};
};

for i = 1 : length(names)
    comp_main(names{i}{1}, mtds, names{i}{2});
end


