function [xyz,xyz_road] = spherical2cartesian(frame)
% converts range, az, el coordinates to x, y, z coordinates for each frame
% input 1 - frame is a cell array of points for each frame in spherical
% coordinates (r, az, el)
% output 1 - xyz is a cell array of points for each frame in x, y, z
% coordinates
% output 2 - xyz_road is a cell array of road/ground points for each frame 
% in x, y, z coordinates
    xyz = cell(0);
    xyz_road = cell(0);
    for i=1:length(frame)
        Q=frame{i};
        %%curate kill all >200m
        Q(Q(:,1)>200,:)=[]; % if range is greater than 200 meters, remove data (all columns)
        r=Q(:,1);
        fee=Q(:,3);
        thet=Q(:,4);
        x=r.*cosd(fee).*cosd(thet);
        y=r.*cosd(fee).*sind(thet);
        z=r.*sind(fee);
        I=Q(:,2);
        I=min(1000,I)/100;
        shot_numbers = Q(:,5);
        points = [x,y,-z,I,shot_numbers];
        xyz{i} = points;
        road_constraints = points(:,1)<200 & points(:,2)>=-3.5 & points(:,2)<=1.5 & points(:,3)>=-2 & points(:,3)<=-1;
        xyz_road{i} = points(road_constraints,:);
    end
end
