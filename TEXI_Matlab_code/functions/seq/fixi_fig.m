function fh = fixi_fig(fig_num, isFEXI, seq, s1, s2, P, opt)

% if nargin < 7
%     lw = 3;
%     fs = 14;
%     ms = 12;
%     YLIM = [.2 1.02];
%     %YLIM = [.1 1.02];
%     opt = {};
% else
%     if isfield(opt,'lw') lw = opt.lw; else lw = 3; end
%     if isfield(opt,'fs') fs = opt.fs; else fs = 14; end
%     if isfield(opt,'ms') ms = opt.ms; else ms = 12; end
%     if isfield(opt,'YLIM') YLIM = opt.YLIM; else YLIM = [.1 1.02]; end
% end

if nargin < 7
    opt.lw = 3;
    opt.fs = 14;
    opt.ms = 12;
else
    if ~isfield(opt,'lw') opt.lw = 3; end
    if ~isfield(opt,'fs') opt.fs = 14; end
    if ~isfield(opt,'ms') opt.ms = 12; end
end

if isFEXI
    series_ind = seq.s_ind;
    bs = seq.b2 * 1e-6;
else
    series_ind = seq.tm_ind;
    bs = seq.b * 1e-6;
end
col = copper(max(series_ind));

if nargin > 4
    RSS = sum((s2(:)-s1(:)).^2);
end
fh = figure(fig_num);
fh.Position(3:4) = [350 400]; %fh.Position(3:4) * 0.7 .* [1 1.2];
fh.Color = 'white';

clf;
hold on
for c = 1:max(series_ind)
    col1 = col(c,:);
    ind = find(series_ind == c);
    b = bs(ind);

    if nargin > 4
        plot(b,s1(ind),'o','LineWidth',opt.lw,'MarkerSize',opt.ms,'color',col1);
        plot(b,s2(ind),'.-','LineWidth',opt.lw,'MarkerSize',2*opt.ms,'color',col1);

    else
        plot(b,s1(ind),'o-','LineWidth',opt.lw,'MarkerSize',opt.ms,'color',col1);
    end
end

set(gca,'YScale', 'log','LineWidth',opt.lw,'Box','off','TickDir','out','TickLength',[.02 .02],'FontSize',opt.fs)
ylabel('signal')
xlabel('b [ms/mm^2]')

if isfield(opt,'XLIM')
    xlim(opt.XLIM)
end
if isfield(opt,'YLIM')
    ylim(opt.YLIM)
end
%pbaspect([1 1 1]);

if isfield(opt, 'showGrid')
    if opt.showGrid
        grid on 
    end
end

if isfield(opt,'title_str')
    title_str{1} = opt.title_str; %sprintf('%s', strrep(opt.seq_name,'_','-'));
end
if exist('P')
    if ~isfield(opt, 'title_relative_fs') opt.title_relative_fs = .8; end

    % check if we have AXR fit
    isAXR = P(2) > 1e-17;

    if isFEXI & isAXR
        title_str{2} = sprintf('Deq = %.2f mm^2/ms, sigma = %.2f, k = %.2f s^{-1}, RSS = %g', P(1)*1e9, P(2), P(3), RSS);
    else
        title_str{2} = sprintf('MD = %.2f mm^2/ms, VD = %.2f mm^4/ms^2, k = %.2f s^{-1}, RSS = %g', P(1)*1e9, P(2)*1e18, P(3), RSS);
    end
end
if exist('title_str')
    th = title(title_str,'FontSize',opt.fs*0.8);
    if isfield(opt, 'title_relative_fs')
        th.FontSize = th.FontSize * opt.title_relative_fs;
    end

    if isfield(opt, 'title_color')
        th.Color = opt.title_color;
    end
end

end