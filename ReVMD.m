function [main_modes, residual_mode, main_modes_hat, residual_hat, ...
    fuzzy_omega, U_matrix, q_scores, keep_scores, ...
    cluster_energy_ratio, u, omega_final] = ...
    ReVMD(signal, maxAlpha, tau, K2, tol, stopc)

    if nargin < 6 || isempty(stopc)
        stopc = 4;
    end

    if nargin < 5 || isempty(tol)
        tol = 1e-6;
    end

    if mod(length(signal), 2) > 0
        signal = signal(1:end-1);
    end

    save_T = length(signal);

    [u, ~, omega_final] = ...
        extractPMFs(signal, maxAlpha, tau, tol, stopc);

    u = ensureModesAsRows(u, save_T);

    omega_final = omega_final(:).';

    K1 = size(u, 1);

    if K1 == 0
        main_modes = [];
        residual_mode = [];
        main_modes_hat = [];
        residual_hat = [];
        fuzzy_omega = [];
        U_matrix = [];
        q_scores = [];
        keep_scores = [];
        cluster_energy_ratio = [];
        return;
    end

    E_i = sum(u.^2, 2).';

    q_scores = E_i / max(E_i + eps);

    if K2 == 1

        fuzzy_omega = ...
            sum(q_scores .* omega_final) / ...
            (sum(q_scores) + eps);

        U_matrix = ones(1, K1);

    else

        idx_init = round(linspace(1, K1, K2));

        centers = omega_final(idx_init).';

        omega_row = omega_final(:).';

        R_matrix = repmat(q_scores, K2, 1);

        max_iter_fcm = 100;

        fcm_tol = 1e-5;

        dist = max(abs(omega_row - centers), eps);

        d_inv = dist.^(-2);

        U_matrix = ...
            d_inv ./ max(sum(d_inv, 1), eps);

        for iter = 1:max_iter_fcm

            Um = U_matrix.^2;

            centers_new = ...
                sum(R_matrix .* Um .* omega_row, 2) ./ ...
                max(sum(R_matrix .* Um, 2), eps);

            dist = ...
                max(abs(omega_row - centers_new), eps);

            d_inv = dist.^(-2);

            U_new = ...
                d_inv ./ max(sum(d_inv, 1), eps);

            if max(abs(U_new(:) - U_matrix(:))) < fcm_tol

                centers = centers_new;

                U_matrix = U_new;

                break;

            end

            centers = centers_new;

            U_matrix = U_new;

        end

        [fuzzy_omega, idx_sort] = ...
            sort(centers(:));

        U_matrix = ...
            U_matrix(idx_sort, :);

    end

    keep_scores = ...
        min(1, q_scores / eps);

    weighted_main = ...
        U_matrix .* keep_scores;

    main_modes = ...
        weighted_main * u;

    residual_mode = ...
        (1 - keep_scores) * u;

    residual_mode = ...
        sum(residual_mode, 1);

    cluster_energy = ...
        sum(main_modes.^2, 2).';

    cluster_energy_ratio = ...
        cluster_energy / ...
        max(sum(cluster_energy), eps);

    main_modes_hat = ...
        fftshift(fft(main_modes, [], 2), 2).';

    residual_hat = ...
        fftshift(fft(residual_mode)).';

end


function modes = ensureModesAsRows(modes, signalLength)

    if isempty(modes)
        return;
    end

    if size(modes, 2) ~= signalLength && ...
            size(modes, 1) == signalLength

        modes = modes.';

    end

    if size(modes, 2) > signalLength

        modes = ...
            modes(:, 1:signalLength);

    end

end


function [u, u_hat_out, omega] = ...
    extractPMFs(signal, maxAlpha, tau, tol, stopc)

    if mod(length(signal), 2) > 0

        signal = signal(1:end-1);

    end

    y = sgolayfilt(signal, 8, 25);

    signoise = signal - y;

    save_T = length(signal);

    T = save_T * 2;

    f_mir = zeros(1, T);

    f_mir_noise = zeros(1, T);

    f_mir(1:T/4) = ...
        signal(save_T/2:-1:1);

    f_mir_noise(1:T/4) = ...
        signoise(save_T/2:-1:1);

    f_mir(T/4+1:3*T/4) = ...
        signal;

    f_mir_noise(T/4+1:3*T/4) = ...
        signoise;

    f_mir(3*T/4+1:T) = ...
        signal(save_T:-1:save_T/2+1);

    f_mir_noise(3*T/4+1:T) = ...
        signoise(save_T:-1:save_T/2+1);

    t = (1:T) / T;

    omega_freqs = ...
        t - 0.5 - 1/T;

    f_hat_full = ...
        fftshift(fft(f_mir));

    f_hat_n_full = ...
        fftshift(fft(f_mir_noise));

    half_idx = ...
        (T/2+1):T;

    omega_half = ...
        omega_freqs(half_idx);

    f_hat_half = ...
        f_hat_full(half_idx);

    f_hat_n_half = ...
        f_hat_n_full(half_idx);

    noisepe = ...
        sum(abs(f_hat_n_half).^2);

    N_max_modes = 50;

    minAlpha = 10;

    alpha = ...
        zeros(1, N_max_modes);

    omega_d_Temp = ...
        zeros(1, N_max_modes);

    u_hat_i_half = ...
        zeros(N_max_modes, length(half_idx));

    polm = ...
        zeros(1, N_max_modes);

    sigerror = ...
        zeros(1, N_max_modes);

    normind = ...
        zeros(1, N_max_modes);

    l = 1;

    SC2 = 0;

    BIC = ...
        zeros(1, N_max_modes);

    N = 300;

    while SC2 ~= 1

        if l == 1

            residual_half = ...
                f_hat_half;

        else

            residual_half = ...
                f_hat_half - ...
                sum(u_hat_i_half(1:l-1, :), 1);

        end

        res_pos = ...
            abs(residual_half);

        [~, idx_pk] = ...
            max(res_pos);

        omega_curr = ...
            omega_half(idx_pk);

        u_hat_curr = ...
            zeros(1, length(half_idx));

        lambda_curr = ...
            zeros(1, length(half_idx));

        Alpha_val = ...
            minAlpha;

        n = 1;

        m = 0;

        bf = 0;

        while Alpha_val < (maxAlpha + 1)

            udiff = ...
                tol + eps;

            Alpha_sq = ...
                Alpha_val^2;

            while udiff > tol && n < N

                omega_diff = ...
                    omega_half - omega_curr;

                omega_diff_sq = ...
                    omega_diff.^2;

                omega_diff_4 = ...
                    omega_diff_sq.^2;

                term1 = ...
                    Alpha_sq * omega_diff_4;

                term2 = ...
                    1 + ...
                    term1 .* ...
                    (1 + ...
                    2 * Alpha_val * omega_diff_sq);

                u_hat_next = ...
                    ( ...
                    residual_half + ...
                    term1 .* u_hat_curr + ...
                    0.5 * lambda_curr ...
                    ) ./ term2;

                mag_sq = ...
                    abs(u_hat_next).^2;

                omega_next = ...
                    (omega_half * mag_sq') / ...
                    (sum(mag_sq) + eps);

                lambda_next = ...
                    lambda_curr + ...
                    tau * ...
                    (residual_half - u_hat_next);

                diff_u = ...
                    u_hat_next - u_hat_curr;

                udiff = ...
                    eps + ...
                    sum(abs(diff_u).^2) / ...
                    (sum(abs(u_hat_curr).^2) + eps);

                u_hat_curr = ...
                    u_hat_next;

                omega_curr = ...
                    omega_next;

                lambda_curr = ...
                    lambda_next;

                n = n + 1;

            end

            if abs(m - log(maxAlpha)) > 1

                m = m + 1;

            else

                m = m + 0.05;

                bf = bf + 1;

            end

            if bf >= 2

                Alpha_val = ...
                    Alpha_val + 1;

            end

            if Alpha_val <= (maxAlpha - 1)

                if bf == 1

                    Alpha_val = ...
                        maxAlpha - 1;

                else

                    Alpha_val = ...
                        exp(m);

                end

                lambda_curr(:) = 0;

                n = 1;

            end

        end

        omega_L_final = ...
            max(omega_curr, 0);

        u_hat_i_half(l, :) = ...
            u_hat_curr;

        omega_d_Temp(l) = ...
            omega_L_final;

        alpha(l) = ...
            Alpha_val;

        if nargin >= 5

            switch stopc

                case 1

                    sigerror(l) = ...
                        sum( ...
                        abs( ...
                        f_hat_half - ...
                        sum(u_hat_i_half(1:l, :), 1) ...
                        ).^2 ...
                        );

                    if l >= 50 || ...
                            sigerror(l) <= noisepe

                        SC2 = 1;

                    end

                case 2

                    sum_u = ...
                        sum(u_hat_i_half(1:l, :), 1);

                    normind(l) = ...
                        sum(abs(sum_u - f_hat_half).^2) / ...
                        (sum(abs(f_hat_half).^2) + eps);

                    if l >= 50 || ...
                            normind(l) < 0.005

                        SC2 = 1;

                    end

                case 3

                    sigerror(l) = ...
                        sum( ...
                        abs( ...
                        f_hat_half - ...
                        sum(u_hat_i_half(1:l, :), 1) ...
                        ).^2 ...
                        );

                    BIC(l) = ...
                        2 * T * log(sigerror(l)) + ...
                        (3 * l) * log(2 * T);

                    if ...
                            ( ...
                            l > 1 && ...
                            BIC(l) > BIC(l-1) ...
                            ) || ...
                            l >= 50

                        SC2 = 1;

                    end

                otherwise

                    weights = ...
                        4 * alpha(l) ./ ...
                        ( ...
                        1 + ...
                        2 * alpha(l) * ...
                        ( ...
                        omega_half - ...
                        omega_d_Temp(l) ...
                        ).^2 ...
                        );

                    temp_val = ...
                        sum( ...
                        weights .* ...
                        u_hat_curr .* ...
                        conj(u_hat_curr) ...
                        );

                    if l < 2

                        polm(l) = ...
                            abs(temp_val);

                        polm_temp = ...
                            polm(l);

                        polm(l) = 1;

                    else

                        polm(l) = ...
                            abs(temp_val) / ...
                            polm_temp;

                    end

                    if ...
                            ( ...
                            l > 1 && ...
                            abs( ...
                            polm(l) - polm(l-1) ...
                            ) < 0.001 ...
                            ) || ...
                            l >= 50

                        SC2 = 1;

                    end

            end

        else

            weights = ...
                4 * alpha(l) ./ ...
                ( ...
                1 + ...
                2 * alpha(l) * ...
                ( ...
                omega_half - ...
                omega_d_Temp(l) ...
                ).^2 ...
                );

            temp_val = ...
                sum( ...
                weights .* ...
                u_hat_curr .* ...
                conj(u_hat_curr) ...
                );

            if l < 2

                polm(l) = ...
                    abs(temp_val);

                polm_temp = ...
                    polm(l);

                polm(l) = 1;

            else

                polm(l) = ...
                    abs(temp_val) / ...
                    polm_temp;

            end

            if ...
                    ( ...
                    l > 1 && ...
                    abs( ...
                    polm(l) - polm(l-1) ...
                    ) < tol ...
                    ) || ...
                    l >= 50

                SC2 = 1;

            end

        end

        l = l + 1;

    end

    l_total = ...
        l - 1;

    omega = ...
        omega_d_Temp(1:l_total);

    u_hat_matrix = ...
        zeros(T, l_total);

    u_hat_matrix((T/2+1):T, :) = ...
        u_hat_i_half(1:l_total, :).';

    u_hat_matrix((T/2+1):-1:2, :) = ...
        conj( ...
        u_hat_matrix((T/2+1):T, :) ...
        );

    u_hat_matrix(1, :) = ...
        conj(u_hat_matrix(end, :));

    u_temp = ...
        real( ...
        ifft( ...
        ifftshift(u_hat_matrix, 1), ...
        [], ...
        1 ...
        ) ...
        ).';

    [omega, indic] = ...
        sort(omega);

    u = ...
        u_temp( ...
        indic, ...
        T/4+1 : 3*T/4 ...
        );

    u_hat_out = ...
        fftshift(fft(u, [], 2), 2).';

end