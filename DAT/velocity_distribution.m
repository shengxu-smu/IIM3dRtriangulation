clear all

load ./Re10/J11/xc.dat
load ./Re10/J11/yc.dat
load ./Re10/J11/zc.dat

nx=size(xc,1);
ny=size(yc,1);
nz=size(zc,1);

% i0=65;
j0=65;
k0=65;
%%
load ./Re10/J11/u.dat
% load ./Re10/J11/v.dat
% load ./Re10/J11/w.dat
% load ./Re10/J11/p.dat
for k=1:nz
    u1(:,:,k)=u(1+(k-1)*ny:k*ny,1:nx);
%     v1(:,:,k)=v(1+(k-1)*ny:k*ny,1:nx);
%     w1(:,:,k)=w(1+(k-1)*ny:k*ny,1:nx);
end
xu1(:) = u1(j0,:,k0);
% yv1(:) = v1(:,i0,k0);
% zw1(:) = w1(i0,j0,:);
%%

load ./Re20/J11/u.dat
% load ./Re20/J11/v.dat
% load ./Re20/J11/w.dat
% load ./Re20/J11/p.dat
for k=1:nz
    u2(:,:,k)=u(1+(k-1)*ny:k*ny,1:nx);
%     v2(:,:,k)=v(1+(k-1)*ny:k*ny,1:nx);
%     w2(:,:,k)=w(1+(k-1)*ny:k*ny,1:nx);
end
xu2(:) = u2(j0,:,k0);
% yv2(:) = v2(:,i0,k0);
% zw2(:) = w2(i0,j0,:);
%%

load ./Re100/J11/u.dat
% load ./Re100/J11/v.dat
% load ./Re100/J11/w.dat
% load ./Re100/J11/p.dat
for k=1:nz
    u3(:,:,k)=u(1+(k-1)*ny:k*ny,1:nx);
%     v3(:,:,k)=v(1+(k-1)*ny:k*ny,1:nx);
%     w3(:,:,k)=w(1+(k-1)*ny:k*ny,1:nx);
end
xu3(:) = u3(j0,:,k0);
% yv3(:) = v3(:,i0,k0);
% zw3(:) = w3(i0,j0,:);
%%
load ./Re200/J11/u.dat
% load ./Re200/J11/v.dat
% load ./Re200/J11/w.dat
% load ./Re200/J11/p.dat
for k=1:nz
    u4(:,:,k)=u(1+(k-1)*ny:k*ny,1:nx);
%     v4(:,:,k)=v(1+(k-1)*ny:k*ny,1:nx);
%     w4(:,:,k)=w(1+(k-1)*ny:k*ny,1:nx);
end
xu4(:) = u4(j0,:,k0);
% yv4(:) = v4(:,i0,k0);
% zw4(:) = w4(i0,j0,:);
%%
figure(1)
C = zeros(1,257)
plot(xc',xu1,'b-.',xc',xu2,'r--',xc',xu3,'g:',xc',xu4,'k-.','LineWidth',1)
legend('Re=10','Re=20','Re=100','Re=200','FontSize', 12)
xlabel('x','FontSize', 12)
ylabel('u','FontSize', 12)
hold on 
plot(xc',C,'k-')
hold off