% This file rearrange date based on their treatments

function data_allocated = Allocate_data(data)
    % determine number of treatments
    num_of_data = length(data);
    names = string([]);

    for cell_idx = 1: num_of_data
        names = [names;data{cell_idx,1}.name(1:find(data{cell_idx,1}.name == '_', 1, 'last')-1)];
    end

    num_of_treatments = length(unique(names));
    unique_names = unique(names);

    num_of_info = length(data{1,1});

    data_allocated = cell(num_of_treatments,1);
    for treat_idx = 1:num_of_treatments
        data_allocated{treat_idx,1}.name = unique_names(treat_idx);
        data_allocated{treat_idx,1}.hetero_radius = [];
        data_allocated{treat_idx,1}.hetero_spacing = [];
        data_allocated{treat_idx,1}.ave_lad_thickness = [];
        data_allocated{treat_idx,1}.lad_thickness_mixed = [];
        data_allocated{treat_idx,1}.lads2total = [];
        data_allocated{treat_idx,1}.g_me = [];
        data_allocated{treat_idx,1}.nucleus_radius = [];
        data_allocated{treat_idx,1}.locs_density = [];

        for cell_idx = 1:num_of_data
            if data{cell_idx,1}.name(1:find(data{cell_idx,1}.name == '_', 1, 'last')-1) == unique_names(treat_idx)
                data_allocated{treat_idx,1}.hetero_radius = ...
                    [data_allocated{treat_idx,1}.hetero_radius;data{cell_idx}.hetero_radius];
                
                data_allocated{treat_idx,1}.hetero_spacing = ...
                    [data_allocated{treat_idx,1}.hetero_spacing;data{cell_idx}.spacing];
               
                data_allocated{treat_idx,1}.ave_lad_thickness = ...
                    [data_allocated{treat_idx,1}.ave_lad_thickness;data{cell_idx}.lad_thickness];

                data_allocated{treat_idx,1}.lad_thickness_mixed = ...
                    [data_allocated{treat_idx,1}.lad_thickness_mixed;data{cell_idx}.lads_seg_thickness];
               
                data_allocated{treat_idx,1}.lads2total = ...
                    [data_allocated{treat_idx,1}.lads2total;data{cell_idx}.lads2total];

                data_allocated{treat_idx,1}.g_me = ...
                    [data_allocated{treat_idx,1}.g_me;data{cell_idx}.g_me];

                data_allocated{treat_idx,1}.nucleus_radius = ...
                    [data_allocated{treat_idx,1}.nucleus_radius;data{cell_idx}.nucleus_radius];

                data_allocated{treat_idx,1}.locs_density = ...
                    [data_allocated{treat_idx,1}.locs_density;data{cell_idx}.locs_density];
            end
        end
    end

