function ax=tidyaxes(ax,fontsize)

if nargin < 2, fontsize = 12; end
if nargin < 1, ax = gca; end
    
set(ax,'FontSize',fontsize);
set(get(ax,'XLabel'),'FontSize',fontsize);
set(get(ax,'YLabel'),'FontSize',fontsize);
set(get(ax,'Title'),'FontSize',fontsize,'FontWeight','bold');

ax.XAxis.LineWidth = 1.5;
ax.YAxis.LineWidth = 1.5;
ax.TickDir         = 'out';

ax.TickLength = [0.02 0.025];

ax.Color           = [1 1 1];
ax.Box             = 'off';
