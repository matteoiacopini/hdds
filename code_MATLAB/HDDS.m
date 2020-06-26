function [p] = HDDS(x,bounds,discrete,Nbins,ramp,inv_scale,radius,legend)
%% HDDS code for plotting an Half-Disk Denisty Strip (HDDS)
if any(isnan(bounds))
   bounds = [min(x),max(x)];
else
   bounds = sort(bounds,'ascend');
end

% compute color shading from probability
if discrete
   breaks = min(min(x),bounds(1)) : 1 : max(max(x),bounds(2));
   ll     = histcounts(x, breaks, 'Normalization', 'probability');
   val    = ll((breaks>=bounds(1)-1) & (breaks<=bounds(2)-1))';
   Nbins  = length(val);
else
   val = ksdensity(x, linspace(bounds(1),bounds(2),Nbins), 'kernel','epanechnikov')';  % val  = normpdf(x,0,1);
end
col = prob_to_col(val,ramp,inv_scale);

% plot HDDS
p=figure('position',[350,500,600,350]);  box on;
stp = 180/Nbins;
for i=1:Nbins
   % NOTE: fills-in the half disk starting from right (upper bound) to left (lower bound),
   % therefore use col(Nbins-i+1,:) instead of col(i,:)
   alpha = linspace((i-1)*stp, i*stp, 100)/180*pi;
   patch([0, cos(alpha)*radius, 0], [0, sin(alpha)*radius, 0], col(Nbins-i+1,:), 'edgealpha',0.01);
end
set(gca,'XTick',[-radius,radius],'XTickLabel',{num2str(bounds(1)),num2str(bounds(2))},'YTickLabel','');
axis([-radius-0.05,radius+0.05, -0.01,radius+0.05]);

if legend
   colormap(ramp);
   c = colorbar('southoutside');
   c.Label.String = 'Density';
end
end