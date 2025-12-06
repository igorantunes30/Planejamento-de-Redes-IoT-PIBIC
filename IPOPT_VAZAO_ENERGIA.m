function test_ipop()
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

pesos = [1, 0; 0.75, 0.25; 0.5, 0.5; 0.25, 0.75; 0.1, 0.9];
Nc_values = Nc_initial:100:Nc_final;
X_all = zeros(length(pesos), length(Nc_values));
Y_all = zeros(length(pesos), length(Nc_values));
Iter_all = zeros(length(pesos), length(Nc_values));

for p_idx = 1:length(pesos)
    pesoR = pesos(p_idx, 1);
    pesoE = pesos(p_idx, 2);
    X = zeros(size(Nc_values));
    Y = zeros(size(Nc_values));

    for idx = 1:length(Nc_values)
        Nc = Nc_values(idx);
        x0 = ones(n_sf,1)/n_sf;

        funcs.objective = @(p) -sum(UTILIDADE_DE_REDE(lambda, p, Nc, b, Toa));
        funcs.gradient = @(p) reshape(real(-UTILIDADE_DE_REDE_GRAD(lambda, p, Nc, b, Toa)), [], 1);
        funcs.constraints = @(p) sum(p) - 1;
        funcs.jacobian = @(p) sparse(1, 1:n_sf, 1, 1, n_sf);
        funcs.jacobianstructure = @() sparse(1, 1:n_sf, 1, 1, n_sf);

        options.lb = zeros(n_sf,1);
        options.ub = ones(n_sf,1);
        options.cl = 0;
        options.cu = 0;
        options.ipopt.hessian_approximation = 'limited-memory';
        options.ipopt.linear_solver = 'mumps';
        options.ipopt.print_level = 0;
        options.ipopt.max_iter = 1000;

        [p1, ~] = ipopt(x0, funcs, options);
        Rmax = UTILIDADE_DE_REDE(lambda, p1, Nc, b, Toa);
        Emin = MODELO_DE_ENERGIA(V, I_id, I_st,I_tx,I_rx,Trx1,Trx2,RD1,RD2, p1, Nc,T, Toa);

        funcs.objective = @(p) sum(MODELO_DE_ENERGIA(V, I_id, I_st,I_tx,I_rx,Trx1,Trx2,RD1,RD2, p, Nc,T, Toa));
        funcs.gradient = @(p) reshape(real(ENERGIA_GRAD(p, Nc, V, I_id, I_st,I_tx,I_rx,Trx1,Trx2,RD1,RD2,T, Toa)), [], 1);
        funcs.jacobian = @(p) sparse(1, 1:n_sf, 1, 1, n_sf);
        funcs.jacobianstructure = @() sparse(1, 1:n_sf, 1, 1, n_sf);

        [p2, ~] = ipopt(x0, funcs, options);
        Emax = MODELO_DE_ENERGIA(V, I_id, I_st,I_tx,I_rx,Trx1,Trx2,RD1,RD2, p2, Nc,T, Toa);
        Rmin = UTILIDADE_DE_REDE(lambda, p2, Nc, b, Toa);

        alfa = Rmax - Rmin;
        beta = Emax - Emin;
        
        funcs.objective = @(p) double(real(-sum((pesoR/alfa) * UTILIDADE_DE_REDE(lambda, p, Nc, b, Toa) - ...
                                                (pesoE/beta) * MODELO_DE_ENERGIA(V, I_id, I_st, I_tx, I_rx, Trx1, Trx2, RD1, RD2, p, Nc, T, Toa))));

        funcs.gradient = @(p) reshape(real(EFF_GRAD(p, lambda, Nc, b, Toa, V, I_id, I_st,I_tx,I_rx,Trx1,Trx2,RD1,RD2,T, pesoR, pesoE, alfa, beta)), 6, 1);
        funcs.jacobian = @(p) sparse(1, 1:n_sf, 1, 1, n_sf);
        funcs.jacobianstructure = @() sparse(1, 1:n_sf, 1, 1, n_sf);

        [p, info] = ipopt(x0, funcs, options);
        Iter_all(p_idx, idx) = info.iter;

        X(idx) = sum(VAZAO(lambda, p, Nc, b, Toa));
        Y(idx) = sum(MODELO_DE_ENERGIA(V, I_id, I_st,I_tx,I_rx,Trx1,Trx2,RD1,RD2, p, Nc,T, Toa));

        fprintf('Peso %d | Nc = %d | Vazão = %.4f | Energia = %.4f | Iter = %d\n', p_idx, Nc, X(idx), Y(idx), Iter_all(p_idx, idx));
    end

    X_all(p_idx, :) = X;
    Y_all(p_idx, :) = Y;
end

save('ipopt_vazao.mat', 'Nc_values', 'X_all', 'pesos');
save('ipopt_energia.mat', 'Nc_values', 'Y_all', 'pesos');
save('ipopt_iteracoes.mat', 'Nc_values', 'Iter_all', 'pesos');

figure;
subplot(3,1,1); hold on;
for k = 1:size(X_all,1), plot(Nc_values, X_all(k,:)); end
grid on; xlabel('Nc'); ylabel('Vazão');
leg = arrayfun(@(k) sprintf('pesoR=%.2f, pesoE=%.2f', pesos(k,1), pesos(k,2)), 1:size(pesos,1), 'UniformOutput', false);
legend(leg, 'Location','best'); hold off;

subplot(3,1,2); hold on;
for k = 1:size(Y_all,1), plot(Nc_values, Y_all(k,:)); end
grid on; xlabel('Nc'); ylabel('Energia');
legend(leg, 'Location','best'); hold off;

subplot(3,1,3); hold on;
for k = 1:size(Iter_all,1), plot(Nc_values, Iter_all(k,:)); end
grid on; xlabel('Nc'); ylabel('Iterações até convergência');
legend(leg, 'Location','best'); hold off;
end

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
R = (log(lambda .* p .* Nc .* b)) - 2 * (TRAFEGO_DE_CARGA(lambda, p, Nc, Toa));
end

function g = UTILIDADE_DE_REDE_GRAD(lambda, p, Nc, b, Toa)
g = zeros(length(p),1);
for i = 1:length(p)
    g(i) = (lambda * Nc * b / p(i)) - (2 * lambda * Nc * Toa(i) / 10000);
end
g = double(real(g));
end

function g = ENERGIA_GRAD(p, Nc, V, I_id, I_st,I_tx,I_rx,Trx1,Trx2,RD1,RD2,T,Toa)
g = zeros(length(p),1);
for i = 1:length(p)
    g(i) = Nc * V * ( ...
        Toa(i) * I_tx + ...
        0.5 * RD1 * I_st + ...
        0.5 * Trx1(i) * I_rx + ...
        0.5 * Toa(i) * I_tx + ...
        0.5 * (RD2 - Trx1(i)) * I_st + ...
        0.5 * (Trx1(i) + Trx2(i)) * I_rx + ...
        0.5 * (T - (Toa(i) + Trx1(i) + RD1)) * I_id + ...
        0.5 * (T - (Toa(i) + RD2 + Trx2(i))) * I_id );
end
g = double(real(g));
end

function g = EFF_GRAD(p, lambda, Nc, b, Toa, V, I_id, I_st,I_tx,I_rx,Trx1,Trx2,RD1,RD2,T, pesoR, pesoE, alfa, beta)
g1 = UTILIDADE_DE_REDE_GRAD(lambda, p, Nc, b, Toa);
g2 = ENERGIA_GRAD(p, Nc, V, I_id, I_st,I_tx,I_rx,Trx1,Trx2,RD1,RD2,T, Toa);
g = (pesoR/alfa) .* g1 - (pesoE/beta) .* g2;
g = g(1:6);
end
