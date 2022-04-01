% Syntax: Baltic_helcom_BEC_analyses.m (script)
%
% Description
% Analyses of helcom data as provided by BEC.
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


%%% <<< EDIT TO COMMIT TO GITHUB >>>

clc; clear;
close all

% Get SSS-satellite data within a given distance each Helcome in-situ
% observation
r = 25; % Radius distance from reference point (in Km)

fg_save = 1; % flag save figures [1]; or not [0];
fg_format = 'png';

path_root = ('/Volumes/Rogue/Data/');
folder_data = ([path_root ...
    'SSS/Baltic/BEC/Validation/indata/Helcom/indata/']);


% Make a log_file to record status of each ARGO-BEC file [2020/01/21]
folder_log = '/Volumes/Rogue/Data/SSS/Baltic/BEC/Validation/indata/';
fn_log = [folder_log 'helcom_MISSING_20200121.txt'];

folder_figs = ['/Volumes/Rogue/scratch/Validation/Helcom/'];
foldercheck_raf(folder_figs)


folder_in = [folder_data ];


fn_in = ([folder_data 'helcom.csv']);

fn_HEAD = ([folder_data 'helcom_header.txt']);

fn_out = ([folder_data 'helcom.mat']);

vars_str = {...
    'Cruise','Station','Type',...
    'DATE','lat','lon','DEPTH',...
    'PRES','TEMP','PSAL'};


vars_str2 = {...
    '[1] Cruise','[2] Station','[3] Type',...
    '[4] DATE','[5] lat','[6] lon','[7] DEPTH',...
    '[9] PRES','[10] TEMP','[11] PSAL'};


vars_vec = [1:7,9:11];

if exist(fn_out,'file') ~= 2
    
    if exist(fn_in,'file') == 2
        data = readtable(fn_in); % read helcome data by BEC (.csv)
        
        % =================
        % READ CSV
        % you will need to change the number
        % of values to match your file %f for numbers and %s for strings.
        ftoread = fn_HEAD;
        fid = fopen(ftoread);
        HEADER = textscan(fid, '%s', 'Delimiter',{','},...
            'TreatAsEmpty',{'[]'},'CommentStyle','//');
        fclose (fid)
        
        HEADER = HEADER{:};
        
        % ===================
        for n = 1:length(vars_vec)
            NVAR = vars_vec(n);
            
            eval([vars_str{n} ' = data{:,' num2str(NVAR) '};']);
        end
        save(fn_out,vars_str{:})
        
    end
    
else
    load(fn_out)
    
end




ind = PRES<=15;

lon (ind) = [];
lat (ind) = [];
PSAL (ind) = [];
TEMP (ind) = [];
Station (ind) = [];
Type (ind) = [];
DATE (ind) = [];
DEPTH (ind) = [];
Cruise (ind) = [];
















