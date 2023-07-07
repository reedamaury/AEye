function [outputArg1,outputArg2] = plot_xyz(xyz, frame_num)
figure; 
scatter3(xyz{frame_num}(:,1),xyz{frame_num}(:,2),xyz{frame_num}(:,3),10,xyz{frame_num}(:,4),'filled'); 
colorbar; xlim([0 80]); zlim([-8 5]); ylim([-10 5]); xlabel("X (m)");ylabel("Y (m)");zlabel("Z (m)"); grid minor;
% ; zlim([-2 5]); view(-90,90)
end

