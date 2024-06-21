% this function calculate the density threshold using all data

function density_threshold = determine_density_threshold(myDir, pct, dss)
    myFiles = dir(fullfile(myDir, '*.txt')); % gets all txt files in struct
    
    density_collection = [];

    for cell_idx = 1:length(myFiles)
        name = extractBefore(myFiles(cell_idx).name,'.txt');
        locs = importdata(myFiles(cell_idx).name);
        locs = locs(:,1:2);
        locs = unique(locs,'rows');
        fprintf('Pre-analysis %s -- Data Size: %d \n', name, length(locs(:,1)))
        % downsampling data -- default dss = 1
        
        down_sample_size = fix(length(locs(:,1))/dss);
        vec = 1:length(locs(:,1));
        rng('default');
        vec = vec(randperm(length(vec)));
        I = vec(1:down_sample_size);
        locs = locs(I,:);

        x = locs(:,1);
        y = locs(:,2);
        
        dt = delaunayTriangulation(x,y);
        [vertices,connections] = voronoiDiagram(dt);
        voronoi_cells = cellfun(@(x) vertices(x,:),connections,'UniformOutput',false);
        voronoi_areas = cellfun(@(x) polyarea(x(:,1),x(:,2)),voronoi_cells,'UniformOutput',false);
        voronoi_areas = vertcat(voronoi_areas{:});
        
        % cal density
        density_temp = 1./voronoi_areas;
        density_collection = [density_collection; density_temp];
    end

    density_threshold = prctile(density_collection,pct);
    fprintf('Density threshold is %d\n', density_threshold)
end

 