% Syntax: COMPUTE_Baltic_SDN_BEC_colocations
%
% Description
%
%
% Inputs
%
% Master script compute colocations of Argo (ref) to satellite grid
% as computed using function Baltic_SDN_BEC_colocations.m. Implementation
% to use BINNED data (new, 20200720)
%
% Options for idata_type
%              [4] Nominal, 
%              [5] Nodal Sampling,
%              [6] Global BEC (v001);
%              [7] CCI+SSS (v01.07);
%              [8] Model (NEMO) data in Baltic
%
% iregion: There are four Baltic+ study regions (see Baltic+ DUM, p. 28):
% [1] Arkona Basin           [ArB] (55°31'30.22"N,  16°16'15.06"E)
% [2] Bothian Sea            [BOS] (61°55'50.19"N,  19°14'46.24"E)
% [3] Gulf of Finland        [GOF] (59°35'55.07"N,  23°07'27.96"E)
% [4] Northern Baltic Proper [NBP] (57°36'10.05"N,  19°48'25.78"E)
% [5] ALL-Baltic             [ALL]
% [6] ALL-OLD                 [ALL-OLD] all Baltic grid points - old code
%
% Version control
% 
% Update: use Baltic_SDN_colocations_v3r0.m [20200715]
% 
% =========================================================================
% Author: rcatany
%
% =========================================================================
clc; clear

%% Script inputs
iyear = 2011:2013;
imonth = 1:12;
%iregion = 1;    % [1] Atlantic; [2] Pacific; [3] ALL-Arctic
% idata_type = 8;
ibasin = 9;      % Basin number: 9 ([9] Baltic, [7] Arctic)
flag_on = 1; % [flag_on= 1] apply flags as (qc = 49 "good value data"); or not [flag_on = 0]
bin_data = 1; % [bin_data = 1] bin REF, or not [bin_data = 0] - go to section 1.4

% % ======================
% define path to Datasets
path_root = ('/Volumes/Rogue/Data/');


%% Get SSS-satellite data within a given distance from each Argo float
r = 0;%25; % Radius distance platform (in Km)
lat_factor = 1.25; % Areas at high-lats are about 25% bigger than at equator

r_str = num2str(r);

r2 = r*lat_factor; % multiply lat_factor to get a more accurate sampling radius

depth_ref = 10; % Reference depth (m), by gral. assumtion 10 m

% get limits of the basin
[xmin,xmax,ymin,ymax, basin_str] = map_lim_raf(ibasin);

iDATA_TYPE = [4,5,7,8];
iREGION = 6%[1:5];




%% set flag_str to notice filename
if flag_on == 1
    flag_str = 'FLAG#49';
        
elseif flag_on == 0
    flag_str = 'FLAG#NO';
    
end

%% set bin_str to notice in filename

if bin_data == 0
    bin_str = 'NOBIN';
elseif bin_data == 1
    bin_str = 'BINNED';
end


%% loop through each dataset (NM, NS, model and CCI+SSS)
for qq = 1:length(iDATA_TYPE)
    idata_type = iDATA_TYPE(qq);
    
    % string of the satellite version
    if idata_type == 1
        version_str = 'plusv3r0';
        
    elseif idata_type == 3
        version_str = 'v2r0';
        
    elseif idata_type == 4
        version_str = 'v1r0_NOM';
        
    elseif idata_type == 5
        version_str = 'v1r0_NS';
        
    elseif idata_type == 6
        version_str = 'v001';
        
    elseif idata_type == 7
        version_str = 'CCI+SSSv1.7';
        
    elseif idata_type == 8
        version_str = 'REANALYSIS_PHY_003_011';
        
    end
    
    
    
    % Loop through each subregion in the Baltic
    for rr = 1:length(iREGION)
        iregion = iREGION(rr);
        
        % Colocatiion fn_out for each Regional study in the Baltic
        if strcmpi(iregion,'ARB') || iregion == 1     % Arkona Basin
            region = 'ARB';
            
        elseif strcmpi(iregion,'BOS')|| iregion == 2  % Bothian Sea
            region = 'BOS';
            
        elseif strcmpi(iregion,'GOF')|| iregion == 3 % Gulf of Finland
            region = 'GOF';
            
        elseif strcmpi(iregion,'NBP')|| iregion == 4 % North Baltic Proper
            region = 'NBP';
            
        elseif strcmpi(iregion,'ALL') || iregion == 5 % all-Baltic region stats
            region = 'ALL';
            
        elseif strcmpi(iregion,'ALL_OLD') || iregion == 6 % all-Baltic region stats old code [including all grid points in Baltic]
            region = 'ALL_OLD';
            
        end
        
        % % =============================
        % Make folder_out for each dataset colocations
        folder_data = [path_root...
            'SSS/' basin_str...
            '/BEC/Validation/indata/SDN/Colocations/'...
            bin_str '/' flag_str '/monthly/'...
            version_str '/' region '/'];
        
        
        folder_out = [path_root...
            'SSS/' basin_str...
            '/BEC/Validation/indata/SDN/Colocations/'...
            bin_str '/' flag_str '/cat_files/'...
            version_str '/' region '/'];
        
        foldercheck_raf(folder_out);
                
        % % =============================
        % % loop through to make colocations
        for yy = 1:length(iyear)
            for mm = 1:length(imonth)
                
                iYEAR = iyear(yy);
                iMONTH = imonth(mm);
                
                [vars_out] =...
                    Baltic_SDN_BEC_colocations_v3r0(iYEAR,iMONTH,iregion,idata_type,flag_on);
                
            end
        end
        
        year_start_str = num2str(iyear(1));
        year_end_str = num2str(iyear(end));
        
        
        % % make concatetanated files (all years in one file)
        % Output variables to save
        ref_vars_out = {...
            'SALT_irange','TEMP_irange','PRES_irange',...
            'lon_irange','lat_irange','time_irange'...
            };
        
        sat_vars_out = {...
            'sss_irange',...
            'sss_error_irange',...
            'lon_sss_irange','lat_sss_irange'};
        
        colocated_vars_out = {...
            'sss_irange_mn',...
            'sss_irange_md',...
            'sss_error_irange_mn'};
        
        % variables to concatenate
        cat_vars = [ref_vars_out sat_vars_out colocated_vars_out];
        
        save_vars_out = [cat_vars {'version_str'}];
        
        ncount = 0; % set counter to zero
        for yy = 1: length(iyear)
            for mm = 1:length(imonth)
                
                ncount = ncount+1;
                
                iYEAR = iyear(yy);
                iMONTH = imonth(mm);
                
                fn_out = [folder_out basin_str '_' region '_SDN_R' r_str '_'...
                    sprintf('%02.0f',depth_ref) 'm_COLOCATIONS_'...
                    year_start_str '_' year_end_str '_'...
                    version_str '.mat'];
                
                fn_out_exist = exist(fn_out,'file');
                
                
                if ncount == 1
                    
                    fn1 = [folder_data basin_str '_' region '_SDN_R' r_str '_'...
                        sprintf('%02.0f',depth_ref) 'm_COLOCATIONS_'...
                        sprintf('%02.0f',iYEAR) sprintf('%02.0f',iMONTH) '_'...
                        version_str '.mat'];
                    
                    if exist(fn1,'file') == 2
                        A = load(fn1);
                        
                    else
                        error([fn1 'does not exist'])
                    end
                    
                    if ~isempty(A.lon_irange)
                    
                    for nn = 1:length(cat_vars)
                        this_param = cat_vars{nn};
                        eval([this_param ' = A.' this_param ';']);
                        
                    end
                    
                    else
                        ncount = 0;
                    end
                    
                elseif ncount > 1
                    fn2 = [folder_data basin_str '_' region '_SDN_R' r_str '_'...
                        sprintf('%02.0f',depth_ref) 'm_COLOCATIONS_'...
                        sprintf('%02.0f',iYEAR) sprintf('%02.0f',iMONTH) '_'...
                        version_str '.mat'];
                    
                    A = load(fn_out);
                    
                    if exist(fn2,'file') == 2
                        B = load(fn2);
                    else
                        error([fn2 'does not exist'])
                    end
                    
                    if ~isempty(B.lon_irange)
                        for nn = 1:length(cat_vars)
                            this_param = cat_vars{nn};
                            eval([this_param ' = cat(2,A.' this_param  ',B.' this_param ');']);
                            
                        end
                    end
                end
    
                % save output only if there is data to be saved
                if ncount ~= 0
                    save(fn_out,save_vars_out{:});
                    
                end
            end
        end
        
        
    end
    
end