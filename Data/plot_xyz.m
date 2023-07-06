function [outputArg1,outputArg2] = plot_xyz(xyz, frame_num)
figure; 
scatter3(xyz{frame_num}(:,1),xyz{frame_num}(:,2),xyz{frame_num}(:,3),10,xyz{frame_num}(:,4),'filled'); 
colorbar; ylim([-5 5]); zlim([-2 5]); view(-90,90); xlabel("X (m)");ylabel("Y (m)");zlabel("Z (m)"); grid minor;
end

