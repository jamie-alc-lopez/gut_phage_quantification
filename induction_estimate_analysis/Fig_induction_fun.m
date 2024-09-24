%% This script generates Fig. 2C

clear;clc

%Parameters
v_div_p = 1e-2; %phage-particle-to-bacteria ratio
B = 10; %burst size
delta = 1; %dilution rate
gamma = delta; %lysis time

%Bound functions
xi_fun{1} = @(psi) (1/B).*v_div_p.*(psi + delta); %From microscopy VLP/bacteria ratio
alt_xi_fun{1} = @(psi) (1/B).*10.*v_div_p.*(psi + delta); %From WGS VLP/bacteria ratio
xi_fun{2} = @(R_minus_1) (R_minus_1).*(gamma + delta)./(B-1 - R_minus_1); %From mgx fold changes

%Combined induction rate bounds
combined_upper = xi_fun{2}(1e-1);
combined_lower = 1e-3;
ylim_vec = [1e-4,1e1];

%Axis vectors
n = 2000;
var_range{1} = logspace(-5,3,n); %psi
var_range{2} = logspace(-3,1,n); %R minus 1
xtick_cell{1} = [1e-5 1e-1 1e3];
xtick_cell{2} = [1e-3 1e-1 1e1];
xlim_cell{1} = [min(var_range{1}),max(var_range{1})];
xlim_cell{2} = [min(var_range{2}),max(var_range{2})];

%Plotting parameters
colors = [220 38 127; 120 94 240]/255; 
combined_color = [129 134 137]/255;
FontSize = 18;
LineWidth1 = 3;
LineWidth2 = 1;
alpha = 0.2;
load('../common_util/font_config.mat');
var_labels = {{'Total phage','adsorption','rate, \psi (day^{-1})'},...
    {'Excess lysogen', 'copy number,','R-1 (dimensionless)'}};
fig_letter = {'C','D'};
figure
pos = get(gcf,'Position');
pos(3) = (2/3)*1.5*pos(3);
pos(4) = 0.85*pos(4);
set(gcf,'Position',pos)
scale_cell = {'log','log'};
ratio_bool = [false false];

%Loop through induction bound functions and plot
for i = 1:length(xi_fun)
    subplot(1,2,i)
    hold on

    if ~ratio_bool(i)
        x = var_range{i};
    else
        x = var_range{i}./xi_fun{i}(var_range{i});
    end

    plot(x,xi_fun{i}(var_range{i}),'-','LineWidth',LineWidth1,...
        'Color',colors(i,:))


    set(gca,'XScale',scale_cell{i},'YScale','log')
    ylim(ylim_vec)
    set(gca,'FontSize',FontSize,'FontName',FontName)
    xlabel(var_labels{i},'Interpreter','tex',...
        'FontSize',LabelFontSize,'FontName',FontName);

    if i == 1
        ylabel({'Induction rate, \xi','(day^{-1})'},'Interpreter','tex',...
            'FontSize',LabelFontSize,'FontName',FontName)
    end

    xvec = [min(x),max(x)];
    set(gca,'TickLabelInterpreter','tex')
    
    xlim(xlim_cell{i})
    xticks([xtick_cell{i}])

    %Make transparent patch around plausible zone
    patch_x = [xlim_cell{i}(1) xlim_cell{i}(2) xlim_cell{i}(2) xlim_cell{i}(1)];
    patch_y = [combined_lower combined_lower combined_upper combined_upper];
    patch(patch_x,patch_y,combined_color,'FaceAlpha',alpha,'EdgeAlpha',alpha);
    text(-0.25,1.02,fig_letter{i},'Interpreter','tex','Units','normalized',...
        'FontSize',PanelFontSize,'FontName',FontName)

    if i == 2
        plot(1e-1,xi_fun{2}(1e-1),'o','MarkerSize',8,'MarkerFaceColor',colors(i,:),...
            'MarkerEdgeColor',colors(i,:));
        text(0.8*1e-1,1.3*xi_fun{2}(1e-1),'Adult','FontSize',GenFontSize,...
            'Interpreter','tex','HorizontalAlignment','right','FontName',FontName)

        %plot(1,xi_fun{2}(0.4),'o','MarkerSize',8,'MarkerFaceColor','w',...
        %    'MarkerEdgeColor',colors(i,:),'LineWidth',2);
        %text(0.8*1,1.3*xi_fun{2}(1),'Infant','FontSize',GenFontSize,...
        %    'Interpreter','tex','HorizontalAlignment','right','FontName',FontName)

    end

end

name_vector = 'plots/Fig_induction_fun.pdf';
exportgraphics(gcf,name_vector,'ContentType','vector')


