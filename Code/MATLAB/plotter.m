

close all;
figure(6)
hold on
%surf(MM, QQ, Z);
q=mesh(xs, ys,h, 100*abs(h-Z)/max(Z(:)));
%q=mesh(flip(xs), flip(ys),h, h);
%q.EdgeColor = '#101010';
%q.EdgeAlpha = 0.2;
%q.EdgeColor = 'none';
f = (100*abs(h-Z)/max(Z(:)));
mean(f(:))
max(f(:))
max(h(:))
ax = gca;
%ax.DataAspectRatio = [diff(get(gca, 'XLim')) diff(get(gca, 'XLim')) 1.3*diff(get(gca, 'ZLim'))];
ax.ZGrid = 'off';
ax.XGrid = 'off';
ax.YGrid = 'off';
ax.Color = 'none';
xlabel('x [mm]');
%ax.XLim = [0,10];
ax.ZLim = [min(h(:)) max(h(:))];
ax.XLim = [0, max(xs)];
ax.YLim = [0, max(ys)];
ylabel('y [mm]');
zlabel('z [mm (approximate)]');
ax.FontName = 'arial';
ax.FontSize = 16;
ax.LineWidth = 1.25;
%set(text,'color','black')
set(gca,'ColorScale','log')
a=colorbar();
%a.Ticks = [0.05,0.1,0.2,0.5,1,2,3];
a.Ticks = [0.03,0.1,0.3,1,3,10,30];
%a.Ticks = [0.01,0.03,0.10,0.30,1.00];
caxis([0.03 30]);
%caxis([0.01 1.5]);
a.Label.String = 'Error [% of max height]';
%a.Label.String = 'z [mm (approximate)]';
view(35,17);
set(gcf,'units','inches','position',[1,1,6.75,5])

set(gca, 'XTick', 0:5:30)
set(gca, 'YTick', 0:5:30)
