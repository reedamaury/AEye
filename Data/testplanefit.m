[frame, xx] = obtainframedata(1716, 1724, "Both_10.10.10.66_%d.csv");
[xyz,xyz_road] = spherical2cartesian(frame);
data = xyz{1}(:,1:3); 
pc = pointCloud(data);
[model, inlierIndices] = pcfitplane(pc, 0.2, [0, 0, 1], 10);  % adjust distanceThreshold as necessary
model.Parameters'
plane1 = select(pc,inlierIndices);
pcshow(plane1)