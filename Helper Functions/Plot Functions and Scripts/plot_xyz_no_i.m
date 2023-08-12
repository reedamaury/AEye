function plot_xyz_no_i(xyz, downrange, color)
    scatter3(xyz(:,1),xyz(:,2),xyz(:,3),15,'filled', "MarkerFaceColor", color); 
    xlim([0 downrange]); zlim([-4 2]); ylim([-10 10]); xlabel("X (m)");ylabel("Y (m)");zlabel("Z (m)"); grid minor;
    % ; zlim([-2 5]); view(-90,90)
end



