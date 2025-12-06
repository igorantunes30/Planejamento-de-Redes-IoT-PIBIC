function NSGAII_AG_VAZAO_ENERGIA


% ---------- PARÂMETROS (iguais ao AG) ----------
n = 80;               % pop_size
max_iter = 200;       % generations
S_min = 7; S_max = 12;
Nc_initial = 500; Nc_final = 4500; step_Nc = 100;

lambda_val = 6; b = 48;
Toa  = [0.1048, 0.1802, 0.3211, 0.5636, 1.0485, 1.9398];
Trx1 = [1.1048, 1.1802, 1.3211, 1.5636, 2.0485, 2.9398];
Trx2 = [2.1048, 2.1802, 2.3211, 2.5636, 3.0485, 3.9398];
V = 3.3; I_tx = 44; I_rx = 10.5; I_st = 1.4; I_id = 0.0015;
RD1 = 1; RD2 = 2; T = 720;

pesos = [1,0; 0.75,0.25; 0.5,0.5; 0.25,0.75; 0.1,0.9];
Nc_values = Nc_initial:step_Nc:Nc_final;

figure;
sgtitle('Otimização com NSGA-II (alinhado ao AG) em Redes LoRa');
subplot(3,1,1); hold on; title('Vazão Total');   xlabel('Número de Nós'); ylabel('Bps');   grid on;
subplot(3,1,2); hold on; title('Energia Total');  xlabel('Número de Nós'); ylabel('Joules'); grid on;
subplot(3,1,3); hold on; title('Convergência');   xlabel('Iterações');    ylabel('Rank-1'); grid on;

for peso_idx = 1:size(pesos,1)
    pesoR = pesos(peso_idx,1); pesoE = pesos(peso_idx,2);

    vazao_total = [];
    energia_total = [];
    convergencia_media = zeros(max_iter,1);

    for Nc = Nc_values
        % Avaliação multiobjetivo: [E, -R] (minimizar ambos)
        eval_vec = @(p) [sum(energia(V,I_id,I_st,I_tx,I_rx,Trx1,Trx2,RD1,RD2,p,Nc,T,Toa)), ...
                         -sum(utilidade(lambda_val,p,Nc,b,Toa))];

        [pop, fit, rank1_hist] = nsga2_simplex(n, S_max-S_min+1, max_iter, eval_vec);

        % Extrair U e E da população final
        U_vals = arrayfun(@(i) sum(utilidade(lambda_val, pop(i,:), Nc, b, Toa)), 1:size(pop,1))';
        E_vals = arrayfun(@(i) sum(energia(V,I_id,I_st,I_tx,I_rx,Trx1,Trx2,RD1,RD2, pop(i,:), Nc, T, Toa)), 1:size(pop,1))';

        % Rmax/Emax/Rmin/Emin estimados a partir da população final
        [~, iR] = max(U_vals);  Rmax = U_vals(iR); Emax = E_vals(iR);
        [~, iE] = min(E_vals);  Emin = E_vals(iE); Rmin = U_vals(iE);

        % Normalização coerente (como no AG)
        alfa = max(Rmax - Rmin, eps);
        beta = max(Emax - Emin, eps);

        % Escalarização ponderada
        eff = (pesoR./alfa).*U_vals - (pesoE./beta).*E_vals;
        [~, ibest] = max(eff);
        best_solution = pop(ibest,:);

        % Métricas no melhor escalarizado
        vazao_total   = [vazao_total,   sum(vazao(lambda_val, best_solution, Nc, b, Toa))];
        energia_total = [energia_total, sum(energia(V,I_id,I_st,I_tx,I_rx,Trx1,Trx2,RD1,RD2, best_solution, Nc, T, Toa))];

        convergencia_media = convergencia_media + rank1_hist(:);
    end

    convergencia_media = convergencia_media / length(Nc_values);

    subplot(3,1,1); plot(Nc_values, vazao_total, '-o', 'DisplayName', sprintf('Peso (%.2f, %.2f)', pesoR, pesoE));
    subplot(3,1,2); plot(Nc_values, energia_total, '-o', 'DisplayName', sprintf('Peso (%.2f, %.2f)', pesoR, pesoE));
    subplot(3,1,3); plot(1:max_iter, convergencia_media, 'DisplayName', sprintf('Peso (%.2f, %.2f)', pesoR, pesoE));

    % salvamento (compatível com execuções múltiplas)
    vazao_nsgaii   = vazao_total;
    energia_nsgaii = energia_total;
    pesos_nsgaii   = pesos;

    if isfile('nsgaii_vazao.mat')
        dados_vazao = load('nsgaii_vazao.mat');
        vazao_nsgaii = [dados_vazao.vazao_nsgaii; vazao_nsgaii];
        pesos_nsgaii = [dados_vazao.pesos_nsgaii; pesos];
    end
    if isfile('nsgaii_energia.mat')
        dados_energia = load('nsgaii_energia.mat');
        energia_nsgaii = [dados_energia.energia_nsgaii; energia_nsgaii];
    end
    save('nsgaii_vazao.mat',   'Nc_values','vazao_nsgaii','pesos_nsgaii');
    save('nsgaii_energia.mat', 'Nc_values','energia_nsgaii','pesos_nsgaii');
end

subplot(3,1,1); legend;
subplot(3,1,2); legend;
subplot(3,1,3); legend;

end

% ======================= NSGA-II NO SIMPLEX (com projeção) =======================
function [pop, fit, rank1_hist] = nsga2_simplex(N, D, max_gen, eval_vec)
pc = 0.9; pm = 1/D; eta_c = 20; eta_m = 20;

% Pop inicial no simplex
pop = rand(N, D); 
pop = pop ./ sum(pop, 2);

fit = zeros(N, 2);
for i = 1:N
    f = eval_vec(pop(i,:)); fit(i,:) = f;
end

rank1_hist = zeros(max_gen,1);

for gen = 1:max_gen
    [fronts, cd] = fast_non_dominated_sort(fit);

    % Seleção por torneio de dominância/crowding
    mating = zeros(N, D);
    for i = 1:N
        a = randi(N); b = randi(N);
        if dominates_sel(fronts, cd, a, b)
            mating(i,:) = pop(a,:);
        else
            mating(i,:) = pop(b,:);
        end
    end

    % Reprodução
    offspring = zeros(N, D); k = 1;
    while k <= N
        p1 = mating(randi(N),:); p2 = mating(randi(N),:);

        if rand < pc
            [c1, c2] = sbx_crossover(p1, p2, eta_c);
        else
            c1 = p1; c2 = p2;
        end

        c1 = poly_mutation(c1, pm, eta_m);
        c2 = poly_mutation(c2, pm, eta_m);

        % PROJEÇÃO NO SIMPLEX (alinhado ao AG)
        c1 = project_to_simplex(c1);
        c2 = project_to_simplex(c2);

        offspring(k,:) = c1;
        if k+1 <= N, offspring(k+1,:) = c2; end
        k = k + 2;
    end

    % Seleção ambiental
    comb = [pop; offspring];
    fitc = zeros(2*N, 2);
    for i = 1:2*N
        f = eval_vec(comb(i,:)); fitc(i,:) = f;
    end

    [fronts_c, cd_c] = fast_non_dominated_sort(fitc);

    new_pop = zeros(N, D); new_fit = zeros(N, 2); filled = 0; fr = 1;
    while filled < N
        idx = find(fronts_c == fr);
        if isempty(idx), break; end
        if filled + numel(idx) <= N
            new_pop(filled+1:filled+numel(idx), :) = comb(idx, :);
            new_fit(filled+1:filled+numel(idx), :) = fitc(idx, :);
            filled = filled + numel(idx);
        else
            [~, ord] = sort(cd_c(idx), 'descend');
            take = N - filled; sel = idx(ord(1:take));
            new_pop(filled+1:N, :) = comb(sel, :);
            new_fit(filled+1:N, :) = fitc(sel, :);
            filled = N;
        end
        fr = fr + 1;
    end

    pop = new_pop; fit = new_fit;
    rank1_hist(gen) = sum(fronts_c == 1);
end
end

function [fronts, cd] = fast_non_dominated_sort(fit)
N = size(fit,1);
S = cell(N,1); n = zeros(N,1); fronts = zeros(N,1);
F{1} = [];
for p = 1:N
    S{p} = []; n(p) = 0;
    for q = 1:N
        if all(fit(p,:) <= fit(q,:)) && any(fit(p,:) < fit(q,:))
            S{p} = [S{p}, q];
        elseif all(fit(q,:) <= fit(p,:)) && any(fit(q,:) < fit(p,:))
            n(p) = n(p) + 1;
        end
    end
    if n(p) == 0
        fronts(p) = 1; F{1} = [F{1}, p];
    end
end
i = 1;
while ~isempty(F{i})
    Q = [];
    for p = F{i}
        for q = S{p}
            n(q) = n(q) - 1;
            if n(q) == 0
                fronts(q) = i + 1; Q = [Q, q];
            end
        end
    end
    i = i + 1; F{i} = Q;
end

cd = zeros(N,1);
maxF = max(fronts);
for i = 1:maxF
    I = find(fronts == i);
    if isempty(I), continue; end
    Ffit = fit(I,:); m = size(Ffit,2); d = zeros(numel(I),1);
    for j = 1:m
        [vals, ord] = sort(Ffit(:,j));
        d(ord(1))  = inf; 
        d(ord(end))= inf;
        vmin = vals(1); vmax = vals(end);
        if vmax > vmin
            for k = 2:numel(I)-1
                d(ord(k)) = d(ord(k)) + (vals(k+1) - vals(k-1)) / (vmax - vmin);
            end
        end
    end
    cd(I) = d;
end
end

function out = dominates_sel(fronts, cd, a, b)
if fronts(a) < fronts(b)
    out = true;
elseif fronts(a) > fronts(b)
    out = false;
else
    out = cd(a) > cd(b);
end
end

function [c1, c2] = sbx_crossover(p1, p2, eta)
u = rand(size(p1)); beta = zeros(size(p1));
beta(u <= 0.5) = (2*u(u <= 0.5)).^(1/(eta+1));
beta(u > 0.5)  = (2 - 2*u(u > 0.5)).^(-1/(eta+1));
c1 = 0.5*((1 + beta).*p1 + (1 - beta).*p2);
c2 = 0.5*((1 - beta).*p1 + (1 + beta).*p2);
c1 = min(max(c1, 0), 1); 
c2 = min(max(c2, 0), 1);
end

function c = poly_mutation(c, pm, eta)
for i = 1:numel(c)
    if rand < pm
        u = rand;
        if u < 0.5
            delta = (2*u)^(1/(eta+1)) - 1;
        else
            delta = 1 - (2*(1-u))^(1/(eta+1));
        end
        c(i) = c(i) + delta;
    end
end
c = min(max(c, 0), 1);
end

% ======================= MÉTRICAS / MODELOS (iguais ao AG) =======================
function Va = vazao(lambda_val, p, Nc, b, Toa)
Va = (lambda_val .* p .* Nc .* b) .* exp(-2 .* trafego(lambda_val, p, Nc, Toa));
end

function En = energia(V, I_id, I_st, I_tx, I_rx, Trx1, Trx2, RD1, RD2, p, Nc, T, Toa)
En = (0.5 .* p .* Nc .* V .* (Toa .* I_tx + RD1 .* I_st + Trx1 .* I_rx) + ...
      0.5 .* p .* Nc .* V .* (Toa .* I_tx + (RD2 - Trx1) .* I_st + (Trx1 + Trx2) .* I_rx)) + ...
     (0.5 .* p .* Nc .* V .* (T - (Toa + Trx1 + RD1)) .* I_id + ...
      0.5 .* p .* Nc .* V .* (T - (Toa + RD2 + Trx2)) .* I_id);
end

function Re = utilidade(lambda_val, p, Nc, b, Toa)
% adiciona eps para evitar -Inf, como no AG
Re = log(lambda_val .* p .* Nc .* b + eps) - 2 .* trafego(lambda_val, p, Nc, Toa);
end

function Gs = trafego(lambda_val, p, Nc, Toa)
Gs = lambda_val .* p .* Nc .* Toa ./ 10000;
end

% ======================= PROJEÇÃO NO SIMPLEX (igual ao AG) =======================
function x = project_to_simplex(y)
y = y(:);
n = numel(y);
u = sort(y, 'descend');
cssv = cumsum(u);
rho = find(u + (1 - cssv) ./ (1:n)' > 0, 1, 'last');
theta = (cssv(rho) - 1) / rho;
x = max(y - theta, 0);
x = (x' / sum(x)); % retorna linha no simplex
end
