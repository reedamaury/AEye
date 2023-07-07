function [outputArg1,outputArg2] = plot_spherical(frame,frame_num)
figure; 
scatter3(frame{frame_num}(:,1),frame{frame_num}(:,4),-1*frame{frame_num}(:,3),10,frame{frame_num}(:,2),'filled'); 
colorbar;  xlabel("X (m)");ylabel("Y (m)");zlabel("Z (m)"); grid minor;
% ylim([-5 5]); zlim([-2 5]); view(-90,90);
end

