%=============================================PARAMETROS====================================================
S_min = 7;
S_max = 12;
Nc_initial = 500;
Nc_final = 4500;
lambda = 6;
b = 48;
Toa = [0.1048, 0.1802, 0.3211, 0.5636, 1.0485, 1.9398]';
Trx1 = [1.1048, 1.1802, 1.3211, 1.5636, 2.0485, 2.9398]';
Trx2 = [2.1048, 2.1802, 2.3211, 2.5636, 3.0485, 3.9398]';
V = 3.3;
I_tx = 44;
I_rx = 10.5;
I_st = 1.4;
I_id = 0.0015;
RD1 = 1;
RD2 = 2;
T = 720;
n_sf = S_max - S_min + 1;

% Pesos para as combinações lineares
pesos = [1, 0; 0.75, 0.25; 0.5, 0.5; 0.25, 0.75; 0.1, 0.9];
labels = {'a=1, b=0', 'a=0.75, b=0.25', 'a=0.5, b=0.5', 'a=0.25, b=0.75', 'a=0.1, b=0.9'};

% Resolução de Nós
% NOTA: Se os outros arquivos (CVX, GA, etc.) usam passo de 1000, 
% altere abaixo para 1000 para garantir compatibilidade nos gráficos.
Nc_values = Nc_initial:100:Nc_final; 

% Inicializa as matrizes de dados
X_all = zeros(length(pesos), length(Nc_values));
Y_all = zeros(length(pesos), length(Nc_values));
EFF_all = zeros(length(pesos), length(Nc_values));
R_all = zeros(length(pesos), length(Nc_values));

% --- CONFIGURAÇÃO GUROBI ---
params = struct();
params.NonConvex = 2;       
params.NumericFocus = 3;    
params.OutputFlag = 0;      
params.FeasibilityTol = 1e-6;
params.OptimalityTol = 1e-6;

fprintf('Iniciando simulação com Gurobi...\n');

% Executa os cálculos para todas as combinações
for p_idx = 1:length(pesos)
    pesoR = pesos(p_idx, 1);
    pesoE = pesos(p_idx, 2);
    
    X = zeros(size(Nc_values)); 
    Y = zeros(size(Nc_values)); 
    EFF_curve = zeros(size(Nc_values));
    R_curve = zeros(size(Nc_values));

    fprintf('Calculando curva %d/%d (a=%.2f, b=%.2f)...\n', p_idx, length(pesos), pesoR, pesoE);

    for idx = 1:length(Nc_values)
        Nc = Nc_values(idx);

        % --- Passo 1: Normalização (Max R e Min E) ---
        [R_best, E_at_Rmax, Rmax_vec, E_at_Rmax_vec] = otimizar_gurobi_Rmax(lambda, Nc, b, Toa, V, I_id, I_st, I_tx, I_rx, Trx1, Trx2, RD1, RD2, T, n_sf, params);
        [E_best, R_at_Emin, E_best_vec, R_at_Emin_vec] = otimizar_gurobi_Emin(lambda, Nc, b, Toa, V, I_id, I_st, I_tx, I_rx, Trx1, Trx2, RD1, RD2, T, n_sf, params);
        
        % Cálculo de alfa e beta (Variações positivas)
        alfa = Rmax_vec - R_at_Emin_vec;
        beta = E_at_Rmax_vec - E_best_vec; 
        
        % Proteção contra valores nulos
        alfa(abs(alfa) < 1e-5) = 1; 
        beta(abs(beta) < 1e-5) = 1; 

        % --- Passo 2: Otimização do Compromisso (EFF) ---
        [p_opt] = otimizar_gurobi_EFF(lambda, Nc, b, Toa, V, I_id, I_st, I_tx, I_rx, Trx1, Trx2, RD1, RD2, T, n_sf, pesoR, pesoE, alfa, beta, params);
        
        % --- Armazena Resultados ---
        p = p_opt;
        X(idx) = sum(VAZAO(lambda, p, Nc, b, Toa));
        Y(idx) = sum(MODELO_DE_ENERGIA(V, I_id, I_st,I_tx,I_rx,Trx1,Trx2,RD1,RD2, p, Nc,T, Toa));
        
        R_vec_opt = UTILIDADE_DE_REDE(lambda, p, Nc, b, Toa);
        E_vec_opt = MODELO_DE_ENERGIA(V, I_id, I_st,I_tx,I_rx,Trx1,Trx2,RD1,RD2, p, Nc , T, Toa);

        R_curve(idx) = sum(R_vec_opt);
        EFF_curve(idx) = sum( (pesoR./alfa) .* R_vec_opt - (pesoE./beta) .* E_vec_opt );
    end
    
    X_all(p_idx, :) = X;
    Y_all(p_idx, :) = Y;
    EFF_all(p_idx, :) = EFF_curve;
    R_all(p_idx, :) = R_curve;
end

%========================= SALVAMENTO DE ARQUIVOS (.MAT) =========================
fprintf('\nSalvando arquivos .mat para leitura no script de comparação...\n');

% 1. Preparar variáveis com os nomes esperados pelo script de leitura
vazao_gurobi   = X_all;
energia_gurobi = Y_all;
pesos_gurobi   = pesos;

% 2. Salvar gurobi_vazao.mat
save('gurobi_vazao.mat', 'Nc_values', 'vazao_gurobi', 'pesos_gurobi');
fprintf('Arquivo salvo: gurobi_vazao.mat\n');

% 3. Salvar gurobi_energia.mat
save('gurobi_energia.mat', 'Nc_values', 'energia_gurobi', 'pesos_gurobi');
fprintf('Arquivo salvo: gurobi_energia.mat\n');

fprintf('Concluído.\n');


%============================ FUNÇÕES OTIMIZAÇÃO GUROBI ============================
function [R_sum, E_sum, R_vec, E_vec, p_sol] = otimizar_gurobi_Rmax(lambda, Nc, b, Toa, V, I_id, I_st, I_tx, I_rx, Trx1, Trx2, RD1, RD2, T, n_sf, params)
    model = struct();
    num_vars = 3 * n_sf; 
    idx_p = 1:n_sf; idx_R_aux = n_sf+1 : 2*n_sf; idx_G_aux = 2*n_sf+1 : 3*n_sf;
    
    epsilon = 1e-4; 
    model.lb = [epsilon * ones(n_sf, 1); -inf(n_sf, 1); -inf(n_sf, 1)]; 
    model.ub = [ones(n_sf, 1); inf(n_sf, 1); inf(n_sf, 1)];

    model.A = sparse(1, idx_p, 1, 1, num_vars);
    model.rhs = 1;
    model.sense = '=';

    model.gencontypes = zeros(1, n_sf);
    model.genconvars = zeros(n_sf, 2);
    for i = 1:n_sf
        model.gencontypes(i) = 0; % LOG
        model.genconvars(i, :) = [idx_p(i), idx_R_aux(i)]; 
    end

    model.modelsense = 'min';
    obj_c = zeros(num_vars, 1);
    obj_c(idx_R_aux) = -1; 
    
    C_G = lambda * Nc / 10000;
    for i = 1:n_sf
        row = size(model.A, 1) + 1;
        model.A(row, idx_G_aux(i)) = 1;
        model.A(row, idx_p(i)) = - C_G * Toa(i);
        model.rhs = [model.rhs; 0];
        model.sense = [model.sense; '='];
        obj_c(idx_G_aux(i)) = 2; 
    end
    
    model.objcon = - sum(log(lambda * Nc * b));
    model.obj = obj_c;
    
    result = gurobi(model, params);
    
    if strcmp(result.status, 'OPTIMAL')
        p_sol = result.x(idx_p);
    else
        p_sol = (1/n_sf) * ones(n_sf, 1); 
    end
    
    R_vec = UTILIDADE_DE_REDE(lambda, p_sol, Nc, b, Toa);
    E_vec = MODELO_DE_ENERGIA(V, I_id, I_st,I_tx,I_rx,Trx1,Trx2,RD1,RD2, p_sol, Nc,T, Toa);
    R_sum = sum(R_vec);
    E_sum = sum(E_vec);
end

function [E_sum, R_sum, E_vec, R_vec, p_sol] = otimizar_gurobi_Emin(lambda, Nc, b, Toa, V, I_id, I_st, I_tx, I_rx, Trx1, Trx2, RD1, RD2, T, n_sf, params)
    model = struct();
    num_vars = n_sf; 
    idx_p = 1:n_sf;
    
    epsilon = 1e-4; 
    model.lb = epsilon * ones(n_sf, 1); 
    model.ub = ones(n_sf, 1);
    
    model.A = sparse(1, idx_p, 1, 1, num_vars);
    model.rhs = 1;
    model.sense = '=';
    model.modelsense = 'min';
    
    const_factor = 0.5 * Nc * V;
    c_txrxst = (Toa .* I_tx + RD1 .* I_st + Trx1 .* I_rx) + (Toa .* I_tx + (RD2 - Trx1) .* I_st + (Trx1 + Trx2) .* I_rx);
    c_idle = (T - (Toa + Trx1 + RD1)).*I_id + (T - (Toa + RD2 + Trx2)).*I_id;
    model.obj = const_factor * (c_txrxst + c_idle);
    
    result = gurobi(model, params);
    
    if strcmp(result.status, 'OPTIMAL')
        p_sol = result.x(idx_p);
    else
        p_sol = (1/n_sf) * ones(n_sf, 1);
    end
    E_vec = MODELO_DE_ENERGIA(V, I_id, I_st,I_tx,I_rx,Trx1,Trx2,RD1,RD2, p_sol, Nc,T, Toa);
    R_vec = UTILIDADE_DE_REDE(lambda, p_sol, Nc, b, Toa);
    E_sum = sum(E_vec);
    R_sum = sum(R_vec);
end

function [p_sol] = otimizar_gurobi_EFF(lambda, Nc, b, Toa, V, I_id, I_st, I_tx, I_rx, Trx1, Trx2, RD1, RD2, T, n_sf, pesoR, pesoE, alfa, beta, params)
    model = struct();
    num_vars = 3 * n_sf;
    idx_p = 1:n_sf; idx_R_aux = n_sf+1 : 2*n_sf; idx_G_aux = 2*n_sf+1 : 3*n_sf;
    
    epsilon = 1e-4;
    model.lb = [epsilon * ones(n_sf, 1); -inf(n_sf, 1); -inf(n_sf, 1)];
    model.ub = [ones(n_sf, 1); inf(n_sf, 1); inf(n_sf, 1)];
    
    model.A = sparse(1, idx_p, 1, 1, num_vars);
    model.rhs = 1;
    model.sense = '=';
    
    model.gencontypes = zeros(1, n_sf);
    model.genconvars = zeros(n_sf, 2);
    for i = 1:n_sf
        model.gencontypes(i) = 0; 
        model.genconvars(i, :) = [idx_p(i), idx_R_aux(i)]; 
    end
    
    model.modelsense = 'min';
    obj_c = zeros(num_vars, 1);
    
    obj_c(idx_R_aux) = - (pesoR ./ alfa);
    
    const_factor = 0.5 * Nc * V;
    c_txrxst = (Toa .* I_tx + RD1 .* I_st + Trx1 .* I_rx) + (Toa .* I_tx + (RD2 - Trx1) .* I_st + (Trx1 + Trx2) .* I_rx);
    c_idle = (T - (Toa + Trx1 + RD1)).*I_id + (T - (Toa + RD2 + Trx2)).*I_id;
    c_E = const_factor * (c_txrxst + c_idle);
    
    obj_c(idx_p) = (pesoE ./ beta) .* c_E;
    
    C_G = lambda * Nc / 10000;
    for i = 1:n_sf
        row = size(model.A, 1) + 1;
        model.A(row, idx_G_aux(i)) = 1;
        model.A(row, idx_p(i)) = - C_G * Toa(i);
        model.rhs = [model.rhs; 0];
        model.sense = [model.sense; '='];
        
        obj_c(idx_G_aux(i)) = (2 * pesoR / alfa(i));
    end
    
    log_constant = log(lambda * Nc * b);
    model.objcon = - sum((pesoR ./ alfa) .* log_constant);
    model.obj = obj_c;
    
    result = gurobi(model, params);
    
    if strcmp(result.status, 'OPTIMAL')
        p_sol = result.x(idx_p);
    else
        p_sol = (1/n_sf) * ones(n_sf, 1);
    end
end

%=============================== FUNÇÕES MATEMÁTICAS ====================================
function G = TRAFEGO_DE_CARGA(lambda, p, Nc, Toa)
    G = lambda .* p .* Nc .* Toa ./10000;
end

function X = VAZAO(lambda, p, Nc, b, Toa)
    X = lambda .* p .* Nc .* b .* exp(-2 .* TRAFEGO_DE_CARGA(lambda, p, Nc, Toa));
end

function E = MODELO_DE_ENERGIA(V, I_id, I_st,I_tx,I_rx,Trx1,Trx2,RD1,RD2, p, Nc,T, Toa)
    E = (0.5 .* p .* Nc .* V .* (Toa .* I_tx + RD1 .* I_st + Trx1 .* I_rx) +...
              0.5 .* p .* Nc .* V .* (Toa .* I_tx + (RD2 - Trx1) .* I_st +...
              (Trx1 + Trx2) .* I_rx)) + (0.5 .* p .* Nc .* V .* (T - (Toa + Trx1 + ...
              RD1)).*I_id +  0.5 .* p .* Nc .* V .* (T - (Toa + RD2 + Trx2)).*I_id);
end

function R = UTILIDADE_DE_REDE(lambda, p, Nc, b, Toa)
    p_safe = p;
    p_safe(p_safe < 1e-9) = 1e-9;
    R = (log(lambda .* p_safe .* Nc .* b)) - 2 * (TRAFEGO_DE_CARGA(lambda, p, Nc, Toa));
end