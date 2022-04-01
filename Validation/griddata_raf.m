function [x_out,y_out] = griddata_raf(ix,iy,qX,qY)
%
% Syntax: [x_out,y_out] = griddata_raf(ix,iy,qX,qY)
%
% Desccription
% Convertert scattered lon/lat data (e.g. Argo float locations) to regular
% grid (e.g. model or satellite grid)
%
% Inputs
% ix: original input (scattered) logitude
% iy: original input (scattered) latitude
% qX: query (regular) longitude
% qY: query (regular) latitude
%
% Outputs
% x_out: output longitude on the regular grid
% y_out: output latitude on the regular grid
%
% =========================================================================

% Check size of input data
if size(ix) ~= size(iy)
    error('input lon and lat must be equal size')
end

if size(qX) ~= size(qY)
    error('input query grid must be of equal size')
end

% Output data has same size input data
[nele] = length(ix); % number elements
x_out = nan(size(ix));
y_out = nan(size(iy));


for nn = 1:nele
    
    DIST_beta = Distance(ix(nn),iy(nn),qX,qY);
    irange_beta = find(DIST_beta == min(DIST_beta(:)));
    
    D = irange_beta;
    
    if length (D) == 1
        x_out(nn) = qX(D);
        y_out(nn) = qY(D);
        
    elseif length(D) > 1
        for tt = 1:length(D)
            x_out(nn) = qX(D(tt));
            y_out(nn) = qY(D(tt));
            
            nn = nn+tt;
            
        end; clear tt
        
    end
    
    
    clear D *_beta
end
