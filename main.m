clear
clc
% 定义一些常量
dx = 1;
dt = 1;
cc = dx/dt;
cs2 = cc^2/3;
n=200/dx;
m=200/dx;
mstep=200000/dt;
uo = 0.1;
rhoo = 1.00;
nu = 1/250;
fprintf('Re = %.2f\n', uo * n / nu);
tau = 3*dt*nu/dx^2 + 0.5;
omega = 1.0/tau;
Nwri = 10;
% 定义变量
rho = ones(n+1,m+1)*rhoo;
f = zeros(9,n + 1,m + 1);
feq = zeros(9,n + 1,m + 1);
u = zeros(n+1,m+1);
v = zeros(n+1,m+1);
velocity = zeros(n+1,m+1);


%define the weight coefficient.
w=[1/9;1/9;1/9;1/9;1/36;1/36;1/36;1/36;4/9];
e = [1,0 ; 0,1 ; -1,0 ; 0,-1 ; 1,1 ; -1,1 ; -1,-1 ; 1,-1 ; 0,0]*cc;


% initialize the distribution function.
for k=1:1:9
    f(k,:,:)=feq_D2Q9(k,rho,u,v,w,e,cs2);
end


%main cycle
tic;
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

    %bottom boundary u=0,v=0 全反弹
    boundary

    %     % top boundary u=u0 Zou-He格式
    %     rhow = f(9,:,end) + f(1,:,end) + f(3,:,end) + 2*(f(2,:,end)+f(5,:,end)+f(6,:,end));
    %     f(4,:,end) = f(2,:,end);
    %     f(7,:,end) = f(5,:,end) - rhow * uo /2/cc + (f(1,:,end) - f(3,:,end))/2;
    %     f(8,:,end) = f(6,:,end) + rhow * uo /2/cc - (f(1,:,end) - f(3,:,end))/2;

    % 计算宏观量
    rho = squeeze(sum(f,1));
    usum=zeros(n+1,m+1);
    vsum=zeros(n+1,m+1);
    for k=1:1:9
        usum = usum + squeeze(f(k,:,:)) * e(k,1);
        vsum = vsum + squeeze(f(k,:,:)) * e(k,2);
    end
    u = usum ./ rho;
    v = vsum ./ rho;


    %   Move at y=m+1， 非平衡态边界条件
    rho(:,m+1) = rho(:,m);
    u(:,m+1) = uo;
    v(:,m+1) = 0;

    if (mod(kk, 1000) == 1)
        delta_u = (u-u1).^2;
        delta_v = (v-u2).^2;
        delta_u = sum(delta_u(:));
        delta_v = sum(delta_v(:));
        sum_u = sum(u(:));
        sum_v = sum(v(:));
        Res_U = sqrt(delta_u+delta_v)/sqrt(sum_u+sum_v);
        fprintf('Iteration %d, Residual (U): %.8f\n', kk, log10(Res_U));
    end

    Feq=ones(9,n+1,2);
    for k=1:1:9
        Feq(k,:,1) = feq_D2Q9(k,rho(:,m+1),u(:,m+1),v(:,m+1),w,e,cs2);
        Feq(k,:,2) = feq_D2Q9(k,rho(:,m),u(:,m),v(:,m),w,e,cs2);
    end
    f(:,:,m+1)=Feq(:,:,1)+f(:,:,m)-Feq(:,:,2);

    %plot_fig
    visua;
end
toc;