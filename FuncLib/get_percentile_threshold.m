function [threshold, storm_data, voronoi_data] = get_percentile_threshold(locs, pct)
    x = locs(:, 1);
    y = locs(:, 2);

    dt = delaunayTriangulation(x, y);
    [vertices, connections] = voronoiDiagram(dt);
    voronoi_cells = cellfun(@(v) vertices(v, :), connections, 'UniformOutput', false);
    voronoi_areas = cellfun(@(cell) polyarea(cell(:, 1), cell(:, 2)), voronoi_cells, 'UniformOutput', false);
    voronoi_areas = vertcat(voronoi_areas{:});

    % Calculate density
    Density = 1 ./ voronoi_areas;
    storm_data = [locs(:, 1:2), Density];
    threshold = prctile(storm_data(:, 3), pct);

    voronoi_data.vertices = vertices;
    voronoi_data.connections = connections;
    voronoi_data.voronoi_areas = voronoi_areas;
end
