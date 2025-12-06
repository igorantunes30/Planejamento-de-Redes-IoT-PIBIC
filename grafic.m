% Carregar os dados dos .MAT
cvx_vazao      = load('cvx_vazao.mat');       % campos: Nc_values, vazao_cvx, pesos_cvx
cvx_energia    = load('cvx_energia.mat');     % campos: Nc_values, energia_cvx, pesos_cvx
fpa_vazao      = load('fpa_vazao.mat');       % campos: Nc_values, vazao_fpa, pesos_fpa
fpa_energia    = load('fpa_energia.mat');     % campos: Nc_values, energia_fpa, pesos_fpa
ag_vazao       = load('ga_vazao.mat');        % campos: Nc_values, vazao_ga,  pesos_ga
ag_energia     = load('ga_energia.mat');      % campos: Nc_values, energia_ga, pesos_ga
nsga_vazao     = load('nsgaii_vazao.mat');    % campos: Nc_values, vazao_nsgaii, pesos_nsgaii
nsga_energia   = load('nsgaii_energia.mat');  % campos: Nc_values, energia_nsgaii, pesos_nsgaii
ip_vazao       = load('ipopt_vazao.mat');     % campos: Nc_values, X_all, pesos
ip_energia     = load('ipopt_energia.mat');   % campos: Nc_values, Y_all, pesos
% --- ADICIONADO GUROBI ---
gurobi_vazao   = load('gurobi_vazao.mat');    % campos: Nc_values, vazao_gurobi, pesos_gurobi
gurobi_energia = load('gurobi_energia.mat');  % campos: Nc_values, energia_gurobi, pesos_gurobi

% Eixo e pesos (usar do CVX como referência)
Nc_values = cvx_vazao.Nc_values;
pesos     = cvx_vazao.pesos_cvx;
num_pesos = size(pesos, 1);

% Seleção de Nc para histogramas (robusto)
Nc_sel_target = [500, 1500, 2500, 3500, 4500];
idx_sel = ismember(Nc_values, Nc_sel_target);
if ~any(idx_sel)
    pick = min(4, numel(Nc_values));
    idx_sel = round(linspace(1, numel(Nc_values), pick));
    Nc_sel = Nc_values(idx_sel);
else
    Nc_sel = Nc_values(idx_sel);
end

% Plot
figure('Name','Comparação Geral (CVX, FPA, GA, NSGA-II, IPOPT, Gurobi)','NumberTitle','off');
t = tiledlayout(num_pesos,4, 'Padding', 'compact', 'TileSpacing', 'compact');
title(t, 'Comparação por Peso: Vazão e Energia + Histogramas');

for i = 1:num_pesos
    label = sprintf('Peso (%.2f, %.2f)', pesos(i,1), pesos(i,2));
    
    % --- Vazão - gráfico de linha ---
    nexttile
    plot(Nc_values, cvx_vazao.vazao_cvx(i,:),       '-o', 'DisplayName', 'CVX'); hold on;
    plot(Nc_values, fpa_vazao.vazao_fpa(i,:),       '-s', 'DisplayName', 'FPA');
    plot(Nc_values, ag_vazao.vazao_ga(i,:),         '-^', 'DisplayName', 'GA');
    plot(Nc_values, nsga_vazao.vazao_nsgaii(i,:),   '-d', 'DisplayName', 'NSGA-II');
    plot(Nc_values, ip_vazao.X_all(i,:),            '-x', 'DisplayName', 'IPOPT');
    plot(Nc_values, gurobi_vazao.vazao_gurobi(i,:), '-*', 'DisplayName', 'Gurobi', 'LineWidth', 1.5); % Gurobi
    
    title(['Vazão - ', label]);
    xlabel('Número de Nós'); ylabel('Vazão (bps)');
    legend('Location','bestoutside'); grid on; hold off;
    
    % --- Vazão - histograma ---
    nexttile
    dados_vazao = [
        cvx_vazao.vazao_cvx(i,idx_sel);
        fpa_vazao.vazao_fpa(i,idx_sel);
        ag_vazao.vazao_ga(i,idx_sel);
        nsga_vazao.vazao_nsgaii(i,idx_sel);
        ip_vazao.X_all(i,idx_sel);
        gurobi_vazao.vazao_gurobi(i,idx_sel); % Gurobi
    ];
    bar(categorical(string(Nc_sel)), dados_vazao'); % categorias em colunas
    title(['Vazão (Hist) - ', label]);
    xlabel('Número de Nós'); ylabel('Vazão (bps)');
    legend({'CVX','FPA','GA','NSGA-II','IPOPT','Gurobi'}, 'Location','bestoutside');
    grid on;
    
    % --- Energia - gráfico de linha ---
    nexttile
    plot(Nc_values, cvx_energia.energia_cvx(i,:),       '-o', 'DisplayName', 'CVX'); hold on;
    plot(Nc_values, fpa_energia.energia_fpa(i,:),       '-s', 'DisplayName', 'FPA');
    plot(Nc_values, ag_energia.energia_ga(i,:),         '-^', 'DisplayName', 'GA');
    plot(Nc_values, nsga_energia.energia_nsgaii(i,:),   '-d', 'DisplayName', 'NSGA-II');
    plot(Nc_values, ip_energia.Y_all(i,:),              '-x', 'DisplayName', 'IPOPT');
    plot(Nc_values, gurobi_energia.energia_gurobi(i,:), '-*', 'DisplayName', 'Gurobi', 'LineWidth', 1.5); % Gurobi
    
    title(['Energia - ', label]);
    xlabel('Número de Nós'); ylabel('Energia (J)');
    legend('Location','bestoutside'); grid on; hold off;
    
    % --- Energia - histograma ---
    nexttile
    dados_energia = [
        cvx_energia.energia_cvx(i,idx_sel);
        fpa_energia.energia_fpa(i,idx_sel);
        ag_energia.energia_ga(i,idx_sel);
        nsga_energia.energia_nsgaii(i,idx_sel);
        ip_energia.Y_all(i,idx_sel);
        gurobi_energia.energia_gurobi(i,idx_sel); % Gurobi
    ];
    bar(categorical(string(Nc_sel)), dados_energia');
    title(['Energia (Hist) - ', label]);
    xlabel('Número de Nós'); ylabel('Energia (J)');
    legend({'CVX','FPA','GA','NSGA-II','IPOPT','Gurobi'}, 'Location','bestoutside');
    grid on;
end

% ============================ EXPORTAÇÃO PARA CSV =============================
outdir = 'csv_export';
if ~exist(outdir,'dir'); mkdir(outdir); end

% Tabelas auxiliares
T_pesos = table((1:num_pesos)', pesos(:,1), pesos(:,2), ...
    'VariableNames', {'peso_idx','a','b'});
writetable(T_pesos, fullfile(outdir,'pesos.csv'));

T_Nc = table(Nc_values(:), 'VariableNames', {'Nc'});
writetable(T_Nc, fullfile(outdir,'Nc_values.csv'));

% ---------- Long format (todas as metodologias) ----------
metodos = {'CVX','FPA','GA','NSGAII','IPOPT','Gurobi'}; % Adicionado Gurobi

% Vazão
rowsV = {};
% Adicionado gurobi_vazao.vazao_gurobi
matsV = {cvx_vazao.vazao_cvx, fpa_vazao.vazao_fpa, ag_vazao.vazao_ga, ...
         nsga_vazao.vazao_nsgaii, ip_vazao.X_all, gurobi_vazao.vazao_gurobi};

for m = 1:numel(metodos)
    M = matsV{m};
    for i = 1:num_pesos
        a = pesos(i,1); b = pesos(i,2);
        for k = 1:numel(Nc_values)
            rowsV(end+1, :) = {metodos{m}, i, a, b, Nc_values(k), M(i,k)}; %#ok<AGROW>
        end
    end
end
T_vazao_all = cell2table(rowsV, 'VariableNames', {'metodo','peso_idx','a','b','Nc','vazao'});
writetable(T_vazao_all, fullfile(outdir,'vazao_all.csv'));

% Energia
rowsE = {};
% Adicionado gurobi_energia.energia_gurobi
matsE = {cvx_energia.energia_cvx, fpa_energia.energia_fpa, ag_energia.energia_ga, ...
         nsga_energia.energia_nsgaii, ip_energia.Y_all, gurobi_energia.energia_gurobi};

for m = 1:numel(metodos)
    M = matsE{m};
    for i = 1:num_pesos
        a = pesos(i,1); b = pesos(i,2);
        for k = 1:numel(Nc_values)
            rowsE(end+1, :) = {metodos{m}, i, a, b, Nc_values(k), M(i,k)}; %#ok<AGROW>
        end
    end
end
T_energia_all = cell2table(rowsE, 'VariableNames', {'metodo','peso_idx','a','b','Nc','energia'});
writetable(T_energia_all, fullfile(outdir,'energia_all.csv'));

% ---------- Subconjunto para histogramas ----------
T_vazao_hist   = T_vazao_all(ismember(T_vazao_all.Nc, Nc_sel), :);
T_energia_hist = T_energia_all(ismember(T_energia_all.Nc, Nc_sel), :);
writetable(T_vazao_hist,   fullfile(outdir,'vazao_hist.csv'));
writetable(T_energia_hist, fullfile(outdir,'energia_hist.csv'));

% ---------- Wide format por método (Nc x pesos) ----------
% CVX
T_cvx_vazao = array2table([Nc_values(:), cvx_vazao.vazao_cvx.'], ...
    'VariableNames', [{'Nc'}, compose('peso_%d', 1:num_pesos)]);
T_cvx_energia = array2table([Nc_values(:), cvx_energia.energia_cvx.'], ...
    'VariableNames', [{'Nc'}, compose('peso_%d', 1:num_pesos)]);
writetable(T_cvx_vazao,   fullfile(outdir,'cvx_vazao_wide.csv'));
writetable(T_cvx_energia, fullfile(outdir,'cvx_energia_wide.csv'));

% FPA
T_fpa_vazao = array2table([Nc_values(:), fpa_vazao.vazao_fpa.'], ...
    'VariableNames', [{'Nc'}, compose('peso_%d', 1:num_pesos)]);
T_fpa_energia = array2table([Nc_values(:), fpa_energia.energia_fpa.'], ...
    'VariableNames', [{'Nc'}, compose('peso_%d', 1:num_pesos)]);
writetable(T_fpa_vazao,   fullfile(outdir,'fpa_vazao_wide.csv'));
writetable(T_fpa_energia, fullfile(outdir,'fpa_energia_wide.csv'));

% GA
T_ga_vazao = array2table([Nc_values(:), ag_vazao.vazao_ga.'], ...
    'VariableNames', [{'Nc'}, compose('peso_%d', 1:num_pesos)]);
T_ga_energia = array2table([Nc_values(:), ag_energia.energia_ga.'], ...
    'VariableNames', [{'Nc'}, compose('peso_%d', 1:num_pesos)]);
writetable(T_ga_vazao,   fullfile(outdir,'ga_vazao_wide.csv'));
writetable(T_ga_energia, fullfile(outdir,'ga_energia_wide.csv'));

% NSGA-II
T_nsga_vazao = array2table([Nc_values(:), nsga_vazao.vazao_nsgaii.'], ...
    'VariableNames', [{'Nc'}, compose('peso_%d', 1:num_pesos)]);
T_nsga_energia = array2table([Nc_values(:), nsga_energia.energia_nsgaii.'], ...
    'VariableNames', [{'Nc'}, compose('peso_%d', 1:num_pesos)]);
writetable(T_nsga_vazao,   fullfile(outdir,'nsgaii_vazao_wide.csv'));
writetable(T_nsga_energia, fullfile(outdir,'nsgaii_energia_wide.csv'));

% IPOPT
T_ip_vazao = array2table([Nc_values(:), ip_vazao.X_all.'], ...
    'VariableNames', [{'Nc'}, compose('peso_%d', 1:num_pesos)]);
T_ip_energia = array2table([Nc_values(:), ip_energia.Y_all.'], ...
    'VariableNames', [{'Nc'}, compose('peso_%d', 1:num_pesos)]);
writetable(T_ip_vazao,   fullfile(outdir,'ipopt_vazao_wide.csv'));
writetable(T_ip_energia, fullfile(outdir,'ipopt_energia_wide.csv'));

% GUROBI (NOVO)
T_gurobi_vazao = array2table([Nc_values(:), gurobi_vazao.vazao_gurobi.'], ...
    'VariableNames', [{'Nc'}, compose('peso_%d', 1:num_pesos)]);
T_gurobi_energia = array2table([Nc_values(:), gurobi_energia.energia_gurobi.'], ...
    'VariableNames', [{'Nc'}, compose('peso_%d', 1:num_pesos)]);
writetable(T_gurobi_vazao,   fullfile(outdir,'gurobi_vazao_wide.csv'));
writetable(T_gurobi_energia, fullfile(outdir,'gurobi_energia_wide.csv'));