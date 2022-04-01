% Syntax: Baltic_ferrybox_BEC_analyses.m (script)
%
% Description
% Analyses of ferrybox data as provided by BEC.
% FerryBox systems measure data from the near-surface layer
% from designated ship routes.
% There are 7-ship route data set (see Baltic+ DUM, p. 15):
% [1] BalticQueen
% [2] FinnMaid
% [3] Romantica
% [4] SiljaSernade
% [5] Transpaper
% [6] Victoria
% 
%
% Input: matlab files created from the original .csv files
%
% Baltic Study [2011-2013]:
% *[1] 
% *[2] 
% *[3]
%
% Data source:
% BEC products for validation
%
%
% current version: v1r0 (2020/01/20)
%
% History
% -
%
%
% =========================================================================
%
% Author: rcatany
%
% =========================================================================

clc; clear;
close all

% Get SSS-satellite data within a given distance from each Argo float
r = 25; % Radius distance from reference point (in Km)


fg_save = 1; % flag save figures [1]; or not [0];
fg_format = 'png';

path_root = ('/Volumes/Rogue/Data/');
folder_data = ([path_root ...
    'SSS/Baltic/BEC/Validation/indata/ferrybox/']);

% ship routes
SHIPS = {'BalticQueen','FinnMaid','Romantica',...
    'SiljaSernade','Transpaper','Victoria'};

% Make a log_file to record status of each ARGO-BEC file [2020/01/21]
folder_log = '/Volumes/Rogue/Data/SSS/Baltic/BEC/Validation/indata/';
fn_log = [folder_log 'ferrybox_MISSING_20200121.txt'];

folder_figs = ['/Volumes/Rogue/scratch/Validation/ferrybox/'];



folder_in = [folder_data SHIPS{1}];

if exist(fn_in,'file') == 2 && exist(fn_out,'file') ~= 2
    
    
    
end



















