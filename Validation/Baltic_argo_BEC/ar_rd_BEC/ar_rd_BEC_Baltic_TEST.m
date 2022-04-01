% ar_rd_BEC_Baltic_TEST.m
%
% Description
% Investigaiton of the information contained in the BEC-ARGO files.
% Platform number appear to be the same in all files
%

clear all
clc

imonth = 6;
iyear = 2013;
ibasin = 9;

DAYS = [1:10];

% B = nan(length(DAYS),10,10)

for nn = 1:length(DAYS)
    iday = DAYS(nn);
    
    [TT] = ar_rd_BEC_Baltic_FUN(iday,imonth,iyear,ibasin);
    
    if ~isempty(TT)
        
        A = [TT.platform TT.ID TT.CYCLE_NUMBER TT.YEAR TT.MONTH TT.DAY]
        
        [i,j] = size(A);
        
        B(nn,1:i,1:j) = A;
        clear TT
    end
    
end

