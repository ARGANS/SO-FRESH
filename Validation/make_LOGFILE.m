function [OutMsg] = make_LOGFILE(fn_out,msg_log)
%
% Syntax: [OutMsg] = make_LOGFILE(folder_out,fn_out)
%
% Description
% Write a log file of the status of files containing data. For example,
% this function is used in Baltic_argo_BEC_analyses, to record the status
% of the ARGO files prepared by BEC for Validation. The status of those
% files are:
% [1] empty: file with no information within it.
% [2] Missing: there is no file for a given day (ie.there should be one)
% [3] Missing variable: files with missing variables.
% 
% Input
% fn_out: filename (whole path) of the logfile.
% msg_log: message to record in the logfile.
% 
% Output
% OutMsg (optional): print output message on the screen
%
% =========================================================================
%
% Author: rcatany
%
% =========================================================================


% file_log = [folder_out '/' 'ARGO_MISSINGt.txt'];
file_log = fn_out;

%msg_log = (['file: ' fn_in ' needs to be downloaded']);
fid = fopen(fullfile(tempdir, file_log), 'a');
if fid == -1
    clc
    warning(msg_log);
    fid = fopen(file_log,'a');
    
end

fid = fopen(file_log,'a');
OutMsg = fprintf(fid,msg_log);