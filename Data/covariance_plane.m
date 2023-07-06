% Load CSV data
[frame, xx] = obtainframedata(1716, 1724, "Both_10.10.10.66_%d.csv");
[xyz,xyz_road] = spherical2cartesian(frame);

abcd = zeros(4, length(frame));
% loop over frames 
for fr = 1:length(frame)
    data = xyz{fr}(:,1:3); 
    pc = pointCloud(data);
    
    % This function uses the M-estimator SAmple Consensus (MSAC) algorithm to find the plane
    % MSAC algorithm is a variant of the RANdom SAmple Consensus (RANSAC) algorithm.
    [model, inlierIndices] = pcfitplane(pc, 0.2, [0, 0, 1], 10);  % adjust distanceThreshold as necessary
    abcd(:,fr) = model.Parameters';

    % Option to plot planes 
    % plane1 = select(pc,inlierIndices);
    % figure(fr)
    % pcshow(plane1)   
end 
param_var = var(abcd,0,2)


