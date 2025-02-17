if (mod(kk, Nwri) == 1)

    u_mag = sqrt(u.^2 + v.^2); % 计算总速度
    subplot(1, 1, 1);
    imagesc(flipud(u_mag'));
    colorbar
    title('Velocity Magnitude');
    axis equal off; drawnow

end