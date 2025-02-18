if (mod(kk, Nwri) == 1)
    u_mag = sqrt(u.^2 + v.^2); % 计算总速度
    subplot(1, 1, 1);
    imagesc(flipud(u_mag'));
    %colorbar;
    title('Velocity Magnitude');
    axis equal off;
    drawnow;

    % 将当前图像帧写入视频
    frame = getframe(gcf); % 获取当前图像
    writeVideo(videoFile, frame); % 将帧写入视频文件
end