pi = 3.14159265359;
%matplotlib inline
close all
R_in = 5; % ~2cm
alpha = 0.4;
l_c = 2.7;

% Capillary length sqrt(gamma/delta pho g)
% for pure water documented at 2.7e-3m


% Plot meniscus
r = linspace(R_in, 2*R_in, 50);
p = linspace(0, 2*pi, 50); % Angles
[R, P] = meshgrid(r, p);
Z = alpha* R_in *(besselk(1,R/l_c)/besselk(1,R_in/l_c)) .* cos(P); % Solution
x = R.*cos(P);
y = R.*sin(P);

A = alpha * R_in;
C = 0;
l1 = 1/l_c;
r1 = 0;
r2 = 0;
r3 = 0;
t1 = 0;
t2 = 0;

mm = linspace(-2*R_in, 2*R_in, 50);
qq = linspace(-2*R_in, 2*R_in, 50);
[MM, QQ] = meshgrid(mm, qq);


% Test of the cartesian form of the eqation with the translations to check
% that it works
%Z2 = cos(r2)*cos(r3)*(A*((MM-t1)*cos(r1)+(QQ-t2)*sin(r1)).*(besselk(1,sqrt((MM-t1).^2+(QQ-t2).^2)*l1)./(sqrt((MM-t1).^2+(QQ-t2).^2)*besselk(1,R_in*l1))))+sin(r3)*MM+sin(r2)*QQ+C;
%Z2 = A*(MM).*(besselk(1,sqrt(MM.^2+QQ.^2)*l1)./(sqrt(MM.^2+QQ.^2)*besselk(1,R_in*l1)));

disp(R_in*sin(alpha));
figure(1)
hold on;
zlim([-2.5 2.5]);
q=surf(x, y, Z);
%q= surf(MM, QQ, Z2);
q.EdgeColor = '#101010';
q.EdgeAlpha = 0.2;
%set(gca,'visible','off')
Ax = gca;
%Ax.ZAxis.Visible = 'off';
Ax.DataAspectRatio = [1 1 0.6];
Ax.ZGrid = 'off';
Ax.XGrid = 'off';
Ax.YGrid = 'off';
Ax.Color = 'none';
xlabel('x [mm]');
ylabel('y [mm]');
zlabel('z [mm]');
Ax.FontName = 'arial';
Ax.FontSize = 10;
Ax.LineWidth = 1.25;

%set(text,'color','black')
a = colorbar();
caxis([-2.5 2.5]);
a.Label.String = 'z [mm]';
view(-9,9);
set(gcf,'units','inches','position',[1,1,8,4.5])