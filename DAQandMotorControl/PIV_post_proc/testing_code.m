%% Testing with Yuanhang's data

clear;

load('test_data_yuanhang.mat');

u = nan([size(compiled_data(1).u),length(compiled_data)]);
v = u;
vort = u;
N = length(compiled_data);

for k = 1:N
    x = flip(compiled_data(k).x,1);
    y = compiled_data(k).y;
    u(:,:,k) = flip(compiled_data(k).u,1);
    v(:,:,k) = flip(compiled_data(k).v,1);
    vort(:,:,k) = flip(compiled_data(k).vort,1);
    figure(1)
    subplot(2,2,1)
    quiver(x,y,u(:,:,k),v(:,:,k),'k');
    axis equal
    pause(0.1)
end

n = size(x,1);
m = size(x,2);

% extracting a foil mask for future use
% Xu = reshape(u,[],N); % reshape snapshots into column vectors
% Xv = reshape(v,[],N);
Xw = reshape(vort,[],N);
% X = [Xu;Xv];
foil_mask = ones(size(Xw));
foil_mask(Xw==0) = 0;

%% Adding artificial noise and artifacts
% NOTE: noise is ADDED to the known data, NOT REPLACED

% % Adding a blob (a.k.a. flume bottom screw reflection)
% x_blob_range = 110:115;%110;
% y_blob_range = 105:110;%135;
% 
% blob_u = nan;%(rand(length(x_blob_range))-1/2).*0.05;
% blob_v = nan;%(rand(length(x_blob_range))-1/2).*0.05;
% 
% % u_noise = u;
% % v_noise = v;
% 
% for k = 1:10%N
%     u_noise(x_blob_range,y_blob_range,k) = u_noise(x_blob_range,y_blob_range,k) + blob_u;
%     v_noise(x_blob_range,y_blob_range,k) = v_noise(x_blob_range,y_blob_range,k) + blob_v;
%     figure(1)
%     subplot(2,2,2)
%     quiver(x,y,u_noise(:,:,k),v_noise(:,:,k),'k');
%     axis equal
%     pause(0.1)
% end
% 
% % Adding an oscillating ellipse (a.k.a. endplate)
% 
% for k = 1:N
%     center_x = 41;
%     center_y = 99 + 10*sin(2*pi*(1/30)*k);
%     th = 0:pi/50:2*pi;
%     x_circle = round(25*cos(th) + center_x);
%     y_circle = round(15*sin(th) + center_y);
%     for j = 1:length(x_circle)
%         x_circ_range = x_circle(j)-1:x_circle(j)+1;
%         y_circ_range = y_circle(j)-1:y_circle(j)+1;
%         % let's just fake the mask
%         u_noise(x_circ_range,y_circ_range,k) = u_noise(x_circ_range,y_circ_range)*nan;% + (rand(3)-1/2).*0.05;
%         v_noise(x_circ_range,y_circ_range,k) = u_noise(x_circ_range,y_circ_range)*nan;% + (rand(3)-1/2).*0.05;
%     end
%     figure(1)
%     subplot(2,2,3)
%     quiver(x,y,u_noise(:,:,k),v_noise(:,:,k),'k');
%     axis equal
%     pause(0.1)
% end
% 
% % Now let's plot the velocity contour
% n = 163;
% m = 183;
% for k = 1:N
%     U_noise = sqrt(u_noise(:,:,k).^2 + v_noise(:,:,k).^2);
%     figure(1)
%     subplot(2,2,4)
%     contourf(x,y,U_noise,'linestyle','none');
%     colormap(brewermap([],'Greys'));
%     axis equal;
% end

% mask = gappy_mask(Xw, 'manual', N, [n,m], 1);
load('i_dont_want_to_make_this_mask_again.mat');
Xnoise = Xw.*mask;

%% Interpolating the marred vectors using a manual mask and GPOD
% This must be used on phase-averaged data

% % Compile data into a spatio-temporal matrix
% Xu_n = reshape(u_noise,[],N); % reshape snapshots into column vectors
% Xv_n = reshape(v_noise,[],N);
% X_n = [Xu_n; Xv_n]; % stack x and y velocity data into matrix
% X_g = X_n;
% 
% % X_g = gappy_mask(X_n, 'manual', N, [n,m], 2); % generates a gappy data matrix with NaN values in place of undesired vectors
% % load('test_gapped_field.mat'); % since we already generated a gappy field
% 
% for k = 1:N
%     figure(2)
%     u_gappy(:,:,k) = reshape(X_g(1:end/2,k),[n,m]);
%     v_gappy(:,:,k) = reshape(X_g(end/2+1:end,k),[n,m]);
%     quiver(x,y,u_gappy(:,:,k),v_gappy(:,:,k),2,'k');
%     title("Gappy flowfield")
%     axis equal
%     pause(0.1)
% end

% Perform Gappy POD to interpolate the missing vectors
[X_i, details] = GPOD(Xnoise);
X_i = foil_mask.*X_i; % multiply by the foil mask

%% Trying to remove the noise using RPCA

[L,S] = RPCA(X_i);
L = foil_mask.*L; % multiply by the foil mask

%% Reshape data
% % From the GPOD interpolated data:
% u_interp = reshape(X_i(1:end/2,:),[n,m,N]);
% v_interp = reshape(X_i(end/2+1:end,:),[n,m,N]);
% 
% % From the gappy data:
% u_gappy = reshape(X_g(1:end/2,:),[n,m,N]);
% v_gappy = reshape(X_g(end/2+1:end,:),[n,m,N]);
% 
% % From the RPCA data:
% 
% u_L = reshape(L(1:end/2,:),[n,m,N]);
% v_L = reshape(L(end/2+1:end,:),[n,m,N]);
vort_i = reshape(X_i,[n,m,N]);
vort_L = reshape(L,[n,m,N]);
vort_g = reshape(Xnoise,[n,m,N]);

%% Plotting

figure()
% quiver(x, y, u(:,:,8), v(:,:,8),2, 'k');
contourf(x,y,vort(:,:,8), 'linestyle','none','levelstep',0.01);
colormap(brewermap([],'-RdBu'));
title('original field')

figure()
% quiver(x, y, u_gappy(:,:,8), v_gappy(:,:,8),2, 'k');
contourf(x,y,vort_g(:,:,8), 'linestyle','none','levelstep',0.01);
colormap(brewermap([],'-RdBu'));
title('gappy field')

figure()
% quiver(x, y, u_interp(:,:,8), v_interp(:,:,8),2, 'k');
contourf(x,y,vort_i(:,:,8), 'linestyle','none','levelstep',0.01);
colormap(brewermap([],'-RdBu'));
title('gappy interpolated field')

figure()
% quiver(x, y, u_L(:,:,8), v_L(:,:,8),2, 'k');
contourf(x,y,vort_L(:,:,8), 'linestyle','none','levelstep',0.01);
colormap(brewermap([],'-RdBu'));
title('RPCA de-noised field')

%% Plotting 2.0

map = brewermap([],'-YlGnBu');

% n = 10;
for n = 1:30
    
    figure(1)
    subplot(1,3,1)
%     U = sqrt(u(:,:,n).^2 + v(:,:,n).^2);
    U = vort(:,:,n);
    contourf(x, y, U, 'linestyle', 'none','levelstep',0.05);
    title('Original Field', 'interpreter', 'latex', 'fontsize', 20);
    colormap(map);
    caxis([-15,15]);
    axis equal
    axis off

    % figure(2)
    subplot(1,3,2)
%     U = sqrt(u_gappy(:,:,n).^2 + v_gappy(:,:,n).^2);
    U = vort_g(:,:,n);
    contourf(x, y, U, 'linestyle', 'none','levelstep',0.05);
    title('Gappy Field', 'interpreter', 'latex', 'fontsize', 20);
    colormap(map);
    caxis([-15,15]);
    axis equal
    axis off

    % figure(3)
    subplot(1,3,3)
%     U = sqrt(u_interp(:,:,n).^2 + v_interp(:,:,n).^2);
    U = vort_i(:,:,n);
    contourf(x, y, U, 'linestyle', 'none','levelstep',0.05);
    title('Interpolated Field', 'interpreter', 'latex', 'fontsize', 20);
    colormap(map);
    caxis([-15,15]);
    axis equal
    axis off

%     % figure(4)
%     subplot(1,4,4)
% %     U = sqrt(u_L(:,:,n).^2 + v_L(:,:,n).^2);
%     U = vort_L(:,:,n);
%     contourf(x, y, U, 'linestyle', 'none','levelstep',0.05);
%     title('RPCA Filtered Field', 'interpreter', 'latex', 'fontsize', 20);
%     colormap(map);
%     caxis([-15,15]);
%     axis equal
%     axis off

    tit = ['frame_',num2str(n),'.png'];
    saveas(gcf,tit)

end

%% Generating an animation

gifname = 'yuanhang_data_GPOD_RPCA_test.gif';

for ii = 1:N
    figname = ['frame_', num2str(ii), '.png'];
    Z = imread(figname);
    [imind,cm] = rgb2ind(Z,256); 
    
    if ii == 1 
        imwrite(imind,cm,gifname,'gif','Loopcount',inf,'DelayTime',1/15); 
    else 
        imwrite(imind,cm,gifname,'gif','WriteMode','append','DelayTime',1/15); 
    end
end

%% Calculating svd

[Ux, Sx, ~] = svd(Xw,'econ');
for i = 1:20
    Sxn(i) = Sx(i,i);
end
Sx = Sxn;
[Ui, Si, ~] = svd(X_i,'econ');
for i = 1:20
    Sin(i) = Si(i,i);
end
Si = Sin;
[UL, SL, ~] = svd(L,'econ');
for i = 1:20
    SLn(i) = SL(i,i);
end
SL = SLn;
     
figure()
subplot(3,2,1)
plot(Sx,'ko','markerfacecolor',[0,0,0])
title('Singular Values Original')
subplot(3,2,2)
plot(Ux(:,1:20))
title('Singular Vectors Original')

subplot(3,2,3)
plot(Si,'ko','markerfacecolor',[0,0,0])
title('Singular Values Interpolated')
subplot(3,2,4)
plot(Ui(:,1:20))
title('Singular Vectors Interpolated')

subplot(3,2,5)
plot(SL,'ko','markerfacecolor',[0,0,0])
title('Singular Values Filtered')
subplot(3,2,6)
plot(UL(:,1:20))
title('Singular Vectors Filtered')
