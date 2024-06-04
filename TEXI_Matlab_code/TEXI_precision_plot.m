
% Plot results from the precision estimation.

clear all
close all

fn_substrate_tag = '07';

results_folder = fullfile('data_TEXI','results',['sim' fn_substrate_tag]);

% scale std according to number of points per protocol
n_data_per_protocol = [42, 24, 12, 8];
std_scale = sqrt(n_data_per_protocol./n_data_per_protocol(1));
% std_scale = [1 1 1 1];


% make table
if (0)
    nk = 1;
    nR = 1;
    nc = 1:4;
    d = 1:2:11;
    Nd = length(d);
    Nc = length(nc);

    for nd = 1:Nd
        %     for nc = 1:Nc
        data_name = sprintf('fitP_SNR100_N1000_inits10_spheres_d%d_VT',d(nd));

        data_path = fullfile(fileparts(fileparts(pwd)), results_folder, data_name);
        display(sprintf('%s',data_path))

        load(data_path,'fitP_data')

        % scale std according to number of points in the protocol
        for m = 1:4
            fitP_data.fitP_std(m,:,:,:,:) = fitP_data.fitP_std(m,:,:,:,:) * std_scale(m);
        end

        Y(nd, :,:) = squeeze(fitP_data.fitP(:,nc,nk,nR,3));
        eY(nd, :,:) = squeeze(fitP_data.fitP_std(:,nc,nk,nR,3));
    end

    % size x protocol x crusher
    sz = size(Y)



    d = 1:2:11;

    c = 0;
    for nc = 1:size(Y, 3)
        for np = 1:size(Y, 2)
            c = c+1;

            for nd = 1:size(Y, 1)

                Yout = Y(nd, np, nc);
                eYout = eY(nd, np, nc);
                YStrings{c,nd} = sprintf('%5.2f ± %3.2f', Yout, eYout);

                if c == 1
                    %                 labels{c} = sprintf('%s %s',protocol_labels{np},crusher_labels{nc});
                    labels{nd} = sprintf('d = %d µm',d(nd));
                end
            end
        end
    end


    T = cell2table(YStrings,'VariableNames', labels)
    writetable(T,'table.xlsx','Sheet',1,'Range','B1')

end

save_fig = 1;

fs = 16;
YLIM = [3 24];
YLIM = [4 18];
% ----
k_true = 10;

d = 3;
nk = 1;
nR = 1;

data_name = sprintf('fitP_SNR100_N1000_inits10_spheres_d%d_VT',d);
data_path = fullfile(fileparts(fileparts(pwd)), results_folder, data_name);
display(sprintf('%s',data_path))

load(data_path,'fitP_data')

protocol = fitP_data.protocol;
fitP = fitP_data.fitP;
fitP_std = fitP_data.fitP_std;

% scale std according to number of points in the protocol
for m = 1:4
    fitP_std(m,:,:,:,:) = fitP_std(m,:,:,:,:) * std_scale(m);
end
qc = fitP_data.qc;
R = fitP_data.R;
[Np, Nc, Nk, NR, ~] = size(fitP);
color_order = copper(Nc);

x = 1:Np;

fh = figure,clf
hold on
plot([0 Np+1],k_true*[1 1],'--','LineWidth',2, 'Color',.5*[1 1 1])

c = 0;
for nc = 1:Nc
    Y = squeeze(fitP(:,nc,nk,nR,3));
    eY = squeeze(fitP_std(:,nc,nk,nR,3));

    errorbar(x+c, Y, eY,'.','LineWidth',3,'MarkerSize',20)
    c = c+.1;
end
grid on
xticks([1:4])
xlim([0.5 Np+1])
ylim(YLIM)
set(gca,'LineWidth',3,'Yscale','Lin','Box','off','TickDir','out','TickLength',[.02 .02],'FontSize',fs)
colororder([0 0 0; color_order]);

if save_fig
    fig_path = fullfile(fileparts(fileparts(pwd)), results_folder);
    mkdir(fig_path)
    fig_path = fullfile(fig_path, [data_name '.png']);

    display(sprintf('saving ... %s',fig_path))
    print(fh, fig_path,'-dpng','-r300');
end


% ----

d = 1:2:11;
Nd = length(d);
color_order = redblue(Nd);

nc = 4;
nk = 1;
nR = 1;
c = 0;
x = 1:Np;

fh = figure,clf
hold on
plot([0 Np+1],k_true*[1 1],'--','LineWidth',2, 'Color',.5*[1 1 1])

data_name = 'fitP_SNR100_N1000_inits10_spheres';

for nd = 1:Nd

    data_path = sprintf('%s_d%d_VT',data_name,d(nd));
    data_path = fullfile(fileparts(fileparts(pwd)), results_folder, data_path);
    display(sprintf('%s',data_path))

    load(data_path,'fitP_data')

    % scale std according to number of points in the protocol
    for m = 1:4
        fitP_data.fitP_std(m,:,:,:,:) = fitP_data.fitP_std(m,:,:,:,:) * std_scale(m);
    end

    Y = squeeze(fitP_data.fitP(:,nc,nk,nR,3));
    eY = squeeze(fitP_data.fitP_std(:,nc,nk,nR,3));

    errorbar(x+c, Y, eY,'.','LineWidth',3,'MarkerSize',20)
    c = c+.1;

end
grid on
xlim([0.5 Np+1])
xticks([1:4])
ylim(YLIM)
% colororder(color_order);
colororder([0 0 0; color_order]);

set(gca,'LineWidth',3,'Yscale','Lin','Box','off','TickDir','out','TickLength',[.02 .02],'FontSize',fs)

if save_fig
    fig_path = fullfile(fileparts(fileparts(pwd)), results_folder);
    mkdir(fig_path)
    fig_path = fullfile(fig_path, [data_name '.png']);

    display(sprintf('saving ... %s',fig_path))
    print(fh, fig_path,'-dpng','-r300');
end






