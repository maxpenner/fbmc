% OQAM-OFDM/FBMC for FeelMaTyc
% Author: Maxim Penner (maxim.penner@ikt.uni-hannover.de)
% Date: 2017-05-01
% svn-repo: https://192.168.1.4/svn/feelmatyc/
% copyright (c) 2017 Maxim Penner

function [samp_rate] = lte(obj, profile_name)

    efficiency_left_guards = 0;

    if strcmp(profile_name, 'bandwidth_1_4')
        
        obj.n_subc      = 128;
        obj.n_guard_l   = 26;
        obj.n_guard_r   = 26;
        obj.deac_DC     = true;
        
        samp_rate = 1.92e6;
        
    elseif strcmp(profile_name, 'bandwidth_2_5')
        
        obj.n_subc      = 256;
        obj.n_guard_l   = 52;
        obj.n_guard_r   = 53;
        obj.deac_DC     = true;
        
        samp_rate = 3.84e6;
        
    elseif strcmp(profile_name, 'bandwidth_5_0')
        
        obj.n_subc      = 512;
        obj.n_guard_l   = 106;
        obj.n_guard_r   = 105;
        obj.deac_DC     = true;
        
        samp_rate = 7.68e6;
        
    elseif strcmp(profile_name, 'bandwidth_10_0')
        
        obj.n_subc      = 1024;
        obj.n_guard_l   = 212;
        obj.n_guard_r   = 211;
        obj.deac_DC     = true;
        
        samp_rate = 15.36e6;
        
    elseif strcmp(profile_name, 'bandwidth_15_0')
        
        obj.n_subc      = 1536;
        obj.n_guard_l   = 318;
        obj.n_guard_r   = 317;
        obj.deac_DC     = true;
        
        samp_rate = 23.04e6;
        
    elseif strcmp(profile_name, 'bandwidth_20_0')
        
        obj.n_subc      = 2048;
        obj.n_guard_l   = 424;
        obj.n_guard_r   = 423;
        obj.deac_DC     = true;
        
        samp_rate = 30.72e6;
        
    elseif strcmp(profile_name, 'bandwidth_1_4_efficient')
        
        obj.n_subc      = 128;
        obj.n_guard_l   = efficiency_left_guards;
        obj.n_guard_r   = efficiency_left_guards;
        obj.deac_DC     = false;
        
        samp_rate = 1.92e6;
        
    elseif strcmp(profile_name, 'bandwidth_2_5_efficient')
        
        obj.n_subc      = 256;
        obj.n_guard_l   = efficiency_left_guards;
        obj.n_guard_r   = efficiency_left_guards;
        obj.deac_DC     = false;
        
        samp_rate = 3.84e6;
        
    elseif strcmp(profile_name, 'bandwidth_5_0_efficient')
        
        obj.n_subc      = 512;
        obj.n_guard_l   = efficiency_left_guards;
        obj.n_guard_r   = efficiency_left_guards;
        obj.deac_DC     = false;
        
        samp_rate = 7.68e6;
        
    elseif strcmp(profile_name, 'bandwidth_10_0_efficient')
        
        obj.n_subc      = 1024;
        obj.n_guard_l   = efficiency_left_guards;
        obj.n_guard_r   = efficiency_left_guards;
        obj.deac_DC     = false;
        
        samp_rate = 15.36e6;
        
    elseif strcmp(profile_name, 'bandwidth_15_0_efficient')
        
        obj.n_subc      = 1536;
        obj.n_guard_l   = efficiency_left_guards;
        obj.n_guard_r   = efficiency_left_guards;
        obj.deac_DC     = false;
        
        samp_rate = 23.04e6;
        
    elseif strcmp(profile_name, 'bandwidth_20_0_efficient')
        
        obj.n_subc      = 2048;
        obj.n_guard_l   = efficiency_left_guards;
        obj.n_guard_r   = efficiency_left_guards;
        obj.deac_DC     = false;
        
        samp_rate = 30.72e6;        
        
    else
        error('No such LTE profile.');
    end

    % reconfigure the fbmc chain
    obj.reconfigure();
end

