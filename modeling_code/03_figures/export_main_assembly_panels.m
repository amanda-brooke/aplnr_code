function outputDir = export_main_assembly_panels(figureRoot)
    %EXPORT_MAIN_ASSEMBLY_PANELS Export consistently sized manuscript axes.
    outputDir = fullfile(figureRoot, 'assembly_panels');
    if ~isfolder(outputDir), mkdir(outputDir); end
    files = dir(fullfile(figureRoot, '*.fig'));
    for i = 1:numel(files)
        source = openfig(fullfile(files(i).folder, files(i).name), 'invisible');
        axesList = findobj(source, 'Type', 'axes');
        positions = zeros(numel(axesList), 2);
        for j = 1:numel(axesList)
            axesList(j).Units = 'normalized';
            positions(j, :) = axesList(j).Position(1:2);
        end
        [~, order] = sortrows([-positions(:, 2), positions(:, 1)]);
        axesList = axesList(order);
        [~, stem] = fileparts(files(i).name);
        for j = 1:numel(axesList)
            fig = figure('Visible', 'off', 'Color', 'w', 'Units', 'inches', 'Position', [1 1 4.3 3.9]);
            ax = copyobj(axesList(j), fig);
            colormap(ax, colormap(axesList(j)));
            set(ax, 'Units', 'inches', 'Position', [0.9 0.9 2.4 2.4], 'FontName', 'Arial', 'FontSize', 16, 'LineWidth', 0.9, 'Color', 'w', 'XColor', 'k', 'YColor', 'k');
            pbaspect(ax, [1 1 1]);
            ax.XLabel.FontSize = 18; ax.YLabel.FontSize = 18;
            ax.Title.FontSize = 18; ax.Title.Color = 'k';
            if contains(stem, 'diffusivity_msd')
                set(ax.YLabel, 'Units', 'normalized', 'Position', [-0.24 0.5 0]);
            end
            if contains(stem, 'density_control') || contains(stem, 'density_shaplnr')
                xlabel(ax, 'Distance from edge (\mum)', 'FontSize', 18);
                ylabel(ax, 'Areal density (1/\mum^2)', 'FontSize', 18);
                yticks(ax, [0 0.8e-5 1.6e-5]);
                ax.YTickLabelMode = 'auto'; ax.YAxis.Exponent = -5;
                if j > 1
                    ylabel(ax, ''); ax.YTick = []; ax.YAxis.Exponent = 0;
                end
                if contains(stem, 'density_control')
                    xlabel(ax, '');
                else
                    title(ax, '');
                end
            end
            isSpeedMap = contains(stem, 'radial_speed_maps');
            if isSpeedMap
                colormap(ax, flipud(magma_colormap(256))); clim(ax, [0 30]);
                isExperimental = contains(stem, 'experimental');
                if ~(isExperimental && j == 1)
                    ylabel(ax, ''); ax.YTick = [];
                end
                if j == 1
                    title(ax, 'Control');
                else
                    title(ax, 'shAPLNR');
                end
                if ~isExperimental && j == 2
                    cb = colorbar(ax); cb.Units = 'inches';
                    cb.Position = [3.4 0.9 0.12 2.4];
                    cb.Ticks = [0 15 30]; cb.FontSize = 16; cb.Color = 'k';
                    cb.Label.String = 'Radial speed (\mum/h)'; cb.Label.FontSize = 18;
                end
            end
            ax.Position = [0.9 0.9 2.4 2.4];
            lines = findobj(ax, 'Type', 'line');
            for k = 1:numel(lines)
                if lines(k).LineWidth >= 2, lines(k).LineWidth = 1.8; end
            end
            if contains(stem, 'inward_velocity') || contains(stem, 'density_control') || contains(stem, 'density_shaplnr')
                named = findobj(ax, 'Type', 'line', '-not', 'DisplayName', '');
                isDensity = contains(stem, 'density_control') || contains(stem, 'density_shaplnr');
                if ~isempty(named) && (~isDensity || j == 1)
                    lg = legend(ax, flipud(named), 'Location', 'northeast', 'Box', 'off', 'FontSize', 14, 'TextColor', 'k');
                    lg.AutoUpdate = 'off';
                    if contains(stem, 'inward_velocity')
                        lg.Location = 'southeast';
                    end
                end
            end
            base = fullfile(outputDir, sprintf('%s_panel%d', stem, j));
            exportgraphics(fig, [base '.png'], 'Resolution', 300, 'BackgroundColor', 'white', 'Padding', 12);
            exportgraphics(fig, [base '.pdf'], 'ContentType', 'vector', 'BackgroundColor', 'white', 'Padding', 12);
            savefig(fig, [base '.fig']); close(fig);
        end
        close(source);
    end
end
