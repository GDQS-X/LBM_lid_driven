clear
clc
% 定义一些常量
dx = 1;
dt = 1;
cc = dx/dt;
cs2 = cc^2/3;
n=250/dx;
m=250/dx;
mstep=200000/dt;
uo = 0.1;
rhoo = 1.00;
nu = 1/4;
fprintf('Re = %.2f\n', uo * n / nu);
tau = 3*dt*nu/dx^2 + 0.5;
omega = 1.0/tau;
Nwri = 10;
% 定义变量
rho = ones(n,m)*rhoo;
u = zeros(n,m);
v = zeros(n,m);

%顶盖驱动
u(1:n-1,m)=uo;

f = zeros(9,n,m);
feq = zeros(9,n,m);

%define the weight coefficient.
w=[1/9;1/9;1/9;1/9;1/36;1/36;1/36;1/36;4/9];
e = [1,0 ; 0,1 ; -1,0 ; 0,-1 ; 1,1 ; -1,1 ; -1,-1 ; 1,-1 ; 0,0]*cc;

%main cycle
tic;
% 设置视频文件
videoFile = VideoWriter('velocity_video.mp4', 'MPEG-4'); % 设置视频文件名和格式
open(videoFile); % 打开视频文件以准备写入
Res_U=1;
for kk=1:1:mstep

    if (mod(kk, 200) == 1)
        u1 = u;
        u2 = v;
    end

    % collision
    for k =1:1:9
        feq(k,:,:)=feq_D2Q9(k,rho,u,v,w,e,cs2);
        f(k,:,:) = omega*feq(k,:,:)+(1-omega)*f(k,:,:);
    end

    %streaming
    evolution

    boundary

    % 计算宏观量
    rho = squeeze(sum(f,1));
    usum=zeros(n,m);
    vsum=zeros(n,m);
    for k=1:1:9
        usum = usum + squeeze(f(k,:,:)) * e(k,1);
        vsum = vsum + squeeze(f(k,:,:)) * e(k,2);
    end
    u = usum ./ rho;
    v = vsum ./ rho;




    %     % top boundary u=u0 Zou-He格式
    %     rhow = f(9,:,end) + f(1,:,end) + f(3,:,end) + 2*(f(2,:,end)+f(5,:,end)+f(6,:,end));
    %     f(4,:,end) = f(2,:,end);
    %     f(7,:,end) = f(5,:,end) - rhow * uo /2/cc + (f(1,:,end) - f(3,:,end))/2;
    %     f(8,:,end) = f(6,:,end) + rhow * uo /2/cc - (f(1,:,end) - f(3,:,end))/2;




    % %   Move at y=m+1， 非平衡态边界条件
    % rho(:,m+1) = rho(:,m);
    % u(:,m+1) = uo;
    % v(:,m+1) = 0;

    if (mod(kk, 1000) == 1)
        delta_u = (u-u1).^2;
        delta_v = (v-u2).^2;
        delta_u = sum(delta_u(:));
        delta_v = sum(delta_v(:));
        sum_u = sum(u(:));
        sum_v = sum(v(:));
        Res_U = sqrt(delta_u+delta_v)/sqrt(sum_u+sum_v);
        fprintf('Iteration %d, Residual (U): %.8f\n', kk, log10(Res_U));

        if log10(Res_U) < -6
            disp('Velocity residual is below threshold, stopping the simulation.');
            break;  % 退出循环
        end
    end

    % Feq=ones(9,n,2);
    % for k=1:1:9
    %     Feq(k,:,1) = feq_D2Q9(k,rho(:,m),u(:,m),v(:,m),w,e,cs2);
    %     Feq(k,:,2) = feq_D2Q9(k,rho(:,m-1),u(:,m-1),v(:,m-1),w,e,cs2);
    % end
    % f(:,:,m)=Feq(:,:,1)+f(:,:,m-1)-Feq(:,:,2);

    %plot_fig
    visua;
end
close(videoFile);
toc;