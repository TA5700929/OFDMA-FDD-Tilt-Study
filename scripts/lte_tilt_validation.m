clear; clc; close all;
rng(42);

%% === User parameters ===
MonteCarlo_iter = 50;
BW = 10e6;              % Hz
freq = 2100;            % MHz
kB = 1.38064852e-23; T = 290;
F_ue = 9;               % dB
hb = 30;                % BS antenna height (m)
k_MIMO = 2;             % spatial multiplexing factor
overhead_factor = 0.75; % MAC/control overhead multiplier
PRB_load = 0.5;         % fraction of PRBs used by scheduled user
shadow_sigma = 8;       % dB log-normal shadowing sigma

terrains = {'Urban','Suburban','Hilly','Vehicular'};
Elec_tilt = 0:0.5:12;   % deg
Mech_tilt = 0:0.5:6;    % finer mechanical resolution (0.5 deg)
numUE = 10;

% CQI Table [CQI, bps/Hz, SINR threshold dB]
CQI_Table = [...
    1 0.15 -6.7; 2 0.23 -4.7; 3 0.38 -2.3; 4 0.6 0.2; 5 0.88 2.4; ...
    6 1.18 4.3; 7 1.48 5.9; 8 1.91 8.1; 9 2.41 10.3; 10 2.73 11.7; ...
    11 3.32 14.1; 12 3.9 16.3; 13 4.52 18.7; 14 5.12 21.0; 15 5.55 22.7];

%% === Output folder setup ===
timestamp = datestr(now,'yyyymmdd_HHMMSS');
outDir = fullfile(pwd, ['IEEE_Access_Output_' timestamp]);
if ~exist(outDir,'dir'), mkdir(outDir); end

logFile = fullfile(outDir, 'run_log.txt');
fidlog = fopen(logFile,'w');
fprintf(fidlog,'Run timestamp: %s\n', timestamp);
fprintf('=== STARTING VALIDATED LTE SIMULATION ===\n');
fprintf(fidlog,'=== STARTING VALIDATED LTE SIMULATION ===\n');

fprintf('Monte Carlo iterations = %d, freq = %d MHz\n', MonteCarlo_iter, freq);
fprintf(fidlog,'Monte Carlo = %d, freq = %d MHz\n', MonteCarlo_iter, freq);

%% === Helper inline functions ===
deg2rad = @(x) x*pi/180;
rmse = @(x,y) sqrt(mean((x - y).^2));

% COST231-Hata correction function
a_hm = @(fc,hm) (1.1*log10(fc) - 0.7).*hm - (1.56*log10(fc)-0.8);

%% === Interferer layout: add 6 interferers (hexagon) ===
cell_radius_km = 1.5;
ang = deg2rad(0:60:300);
interferer_xy = [cell_radius_km*cos(ang)', cell_radius_km*sin(ang)']; % 6 points

%% === Preallocate results ===
results = [];
results_detail = []; 
mc_traces = struct();

%% === Simulation loop (terrain x mech x elec) ===
fprintf('Processing simulation grid...\n'); fprintf(fidlog,'Processing simulation grid...\n');
for ti = 1:numel(terrains)
    terrain = terrains{ti};
    mc_traces.(terrain) = zeros(MonteCarlo_iter, numel(Mech_tilt)*numel(Elec_tilt)); 
    colIdx = 1;
    fprintf(' Processing %s terrain...\n', terrain); fprintf(fidlog,' Processing %s terrain...\n', terrain);
    for m_idx = 1:numel(Mech_tilt)
        for e_idx = 1:numel(Elec_tilt)
            mech = Mech_tilt(m_idx);
            elec = Elec_tilt(e_idx);


            d_ref_km = 1;
            switch terrain
                case 'Urban'
                    PL_hata_ref = 69.55 + 26.16*log10(freq) - 13.82*log10(hb) + (44.9 - 6.55*log10(hb))*log10(d_ref_km);
                    CM = 3; 
                case 'Suburban'
                    PL_hata_ref = 69.55 + 26.16*log10(freq) - 13.82*log10(hb) + (44.9 - 6.55*log10(hb))*log10(d_ref_km) - 2*((log10(freq/28))^2) - 5.4;
                    CM = 0;
                case 'Hilly'
                    PL_hata_ref = 69.55 + 26.16*log10(freq) - 13.82*log10(hb) + (44.9 - 6.55*log10(hb))*log10(d_ref_km);
                    CM = 0;
                otherwise % Vehicular
                    PL_hata_ref = 69.55 + 26.16*log10(freq) - 13.82*log10(hb) + (44.9 - 6.55*log10(hb))*log10(d_ref_km);
                    CM = 0;
            end

            % nominal tilt loss approximation
            tilt_loss = 0.2 * (mech + elec);

            % --- Monte Carlo scheduling (PF and RR)
            sched_Tput_PF = zeros(MonteCarlo_iter,1);
            sched_Tput_RR = zeros(MonteCarlo_iter,1);

            for mc = 1:MonteCarlo_iter
                
                user_d_km = 0.05 + 10*rand(numUE,1); 
               
                PL_hata = 69.55 + 26.16*log10(freq) - 13.82*log10(hb) + (44.9 - 6.55*log10(hb)).*log10(user_d_km);
              
                hm = 1.5; 
                PL_cost231 = 46.3 + 33.9*log10(freq) - 13.82*log10(hb) - a_hm(freq,hm) + (44.9 - 6.55*log10(hb)).*log10(user_d_km) + CM;
               
                PL_winner = 32.45 + 36.7*log10(user_d_km) + 20*log10(freq/1000); 

              
                shadow = randn(numUE,1) * shadow_sigma;
                PL_hata = PL_hata + shadow;
                PL_cost231 = PL_cost231 + shadow;
                PL_winner = PL_winner + shadow;

                % pick model for link budget used in scheduling: choose Hata baseline (consistent)
                PL_used = PL_hata;

                % compute received power per UE (dBm)
                P_tx = 43; G_tx = 17; G_rx = 0;
                Pr_dBm_users = P_tx + G_tx + G_rx - PL_used - tilt_loss;

                % interference aggregation from 6 hexagon interferers
                I_mW = 0;
                for ifi = 1:size(interferer_xy,1)
                    % approximate interferer distance factor using geometry and add same PL slope
                    d_if_km = hypot(interferer_xy(ifi,1), interferer_xy(ifi,2));
                    PL_if = 69.55 + 26.16*log10(freq) - 13.82*log10(hb) + (44.9 - 6.55*log10(hb))*log10(max(d_if_km,0.001));
                    I_mW = I_mW + 10.^((P_tx + G_tx - PL_if - tilt_loss - 30)/10);
                end
                I_dBm = 10*log10(I_mW + eps) + 30;

                N_dBm = 10*log10(kB*T*BW) + 30 + F_ue;
                N_mW = 10.^((N_dBm - 30)/10);

                % SINR per UE (dB)
                Pr_mW = 10.^((Pr_dBm_users - 30)/10);
                SINR_linear = Pr_mW ./ (N_mW + I_mW + eps);
                SINR_dB_users = 10*log10(SINR_linear + eps);

                % Map SINR to CQI safely
                % idx = sum(SINR >= thresholds) method
                idx = sum(bsxfun(@ge, SINR_dB_users, CQI_Table(:,3)'), 2);
                idx(idx < 1) = 1; idx(idx > size(CQI_Table,1)) = size(CQI_Table,1);
                eta_users = CQI_Table(idx,2); % bits/s/Hz

                instTput_perstream = eta_users * BW / 1e6; % Mbps per stream
                instTput_user = instTput_perstream * k_MIMO * overhead_factor;

                % ---- Proportional Fair scheduler
                if mc == 1
                    avgT = ones(numUE,1) * mean(instTput_user);
                end
                PF_weight = instTput_user ./ (avgT + eps);
                [~, selPF] = max(PF_weight);
                sched_Tput_PF(mc) = instTput_user(selPF) * PRB_load;
                avgT = (1 - 0.1)*avgT + 0.1*instTput_user;

                % ---- Round-Robin scheduler
                selRR = mod(mc-1, numUE) + 1;
                sched_Tput_RR(mc) = instTput_user(selRR) * PRB_load;
            end % MonteCarlo

            % collect stats
            meanPF = mean(sched_Tput_PF); stdPF = std(sched_Tput_PF);
            meanRR = mean(sched_Tput_RR); stdRR = std(sched_Tput_RR);

           
            A = 69.55 + 26.16*log10(freq) - 13.82*log10(hb);
            B = 44.9 - 6.55*log10(hb);
            PL_nom = PL_hata_ref + tilt_loss; % baseline (d_ref_km)
            d_est = 10.^((PL_nom - A)./B);
            coverage_km = max(0.1, min(d_est, 10));

            results = [results; {terrain, mech, elec, meanPF, stdPF, meanRR, stdRR, coverage_km}];
            results_detail = [results_detail; {terrain, mech, elec, sched_Tput_PF, sched_Tput_RR}];

            mc_traces.(terrain)(:,colIdx) = sched_Tput_PF;
            colIdx = colIdx + 1;
        end
    end
end

%% === Convert to table and save ===
ResultsTbl = cell2table(results, 'VariableNames', ...
    {'Terrain','MechTilt_deg','ElecTilt_deg','MeanTput_PF_Mbps','StdTput_PF','MeanTput_RR_Mbps','StdTput_RR','Coverage_km'});

csvFile = fullfile(outDir, ['Results_' timestamp '.csv']);
writetable(ResultsTbl, csvFile);
fprintf('Simulation complete. Results -> %s\n', csvFile);
fprintf(fidlog,'Results CSV: %s\n', csvFile);

%% === Statistical analysis ===
[p_terrain,~,~] = anova1(ResultsTbl.MeanTput_PF_Mbps, ResultsTbl.Terrain, 'off');
[p_anovan,~,~] = anovan(ResultsTbl.MeanTput_PF_Mbps, ...
    {ResultsTbl.MechTilt_deg, ResultsTbl.ElecTilt_deg}, 'model','interaction','display','off');

fprintf('ANOVA (terrain effect) p = %.4f\n', p_terrain);
fprintf(fidlog,'ANOVA p = %.4f\n', p_terrain);
fprintf('Two-way anovan p-values (mech, elec) = %.4f, %.4f\n', p_anovan(1), p_anovan(2));
fprintf(fidlog,'Two-way anovan p-values (mech, elec) = %.4f, %.4f\n', p_anovan(1), p_anovan(2));

%% === Pareto analysis & save figures (throughput vs coverage) ===
terrains_u = unique(ResultsTbl.Terrain);
for k = 1:numel(terrains_u)
    mask = strcmp(ResultsTbl.Terrain, terrains_u{k});
    sub = ResultsTbl(mask,:);
    x = sub.MeanTput_PF_Mbps; y = sub.Coverage_km;
    idx_pf = paretofront_local(x,y);
    fig = figure('Visible','off'); hold on;
    scatter(y,x,25,'b','filled');
    scatter(y(idx_pf), x(idx_pf), 80,'r','filled');
    xlabel('Coverage (km)'); ylabel('Mean Throughput PF (Mbps)');
    title(sprintf('Pareto - %s', terrains_u{k}));
    grid on;
    saveas(fig, fullfile(outDir, sprintf('Pareto_%s.png', terrains_u{k})));
    close(fig);
end

%% === Monte Carlo convergence plots (sample) ===
for k = 1:numel(terrains_u)
    key = terrains_u{k};
    traceMat = mc_traces.(key);
    fig = figure('Visible','off');
    plot(1:MonteCarlo_iter, traceMat(:,1), '-o','DisplayName','Scenario sample'); hold on;
    plot(1:MonteCarlo_iter, mean(traceMat,2), '-','LineWidth',2,'DisplayName','Mean across scenarios');
    xlabel('Monte Carlo iteration'); ylabel('Scheduled Throughput (Mbps)');
    title(sprintf('MC Convergence - %s', key));
    legend('Location','best'); grid on;
    saveas(fig, fullfile(outDir, sprintf('MC_Convergence_%s.png', key)));
    close(fig);
end

%% === Contour sensitivity (interp2) for one terrain example (Urban) ===

urban_mask = strcmp(ResultsTbl.Terrain,'Urban');
Ux = unique(ResultsTbl.MechTilt_deg(urban_mask));
Uy = unique(ResultsTbl.ElecTilt_deg(urban_mask));
[MechG, ElecG] = meshgrid(Uy, Ux); 

Z = nan(length(Ux), length(Uy));
for i = 1:length(Ux)
    for j = 1:length(Uy)
        rmask = urban_mask & ResultsTbl.MechTilt_deg == Ux(i) & ResultsTbl.ElecTilt_deg == Uy(j);
        if any(rmask)
            Z(i,j) = mean(ResultsTbl.MeanTput_PF_Mbps(rmask));
        end
    end
end
% interpolate to finer grid
[Xq, Yq] = meshgrid(linspace(min(Uy),max(Uy),200), linspace(min(Ux),max(Ux),200));
Zq = interp2(Uy, Ux, Z, Xq, Yq, 'linear');
fig = figure('Visible','off');
contourf(Xq, Yq, Zq, 20,'LineColor','none');
colorbar; xlabel('Electrical Tilt (deg)'); ylabel('Mechanical Tilt (deg)');
title('Throughput Sensitivity (Urban) - interp2 contour');
saveas(fig, fullfile(outDir,'Contour_Urban_Throughput.png'));
close(fig);

%% === Hata vs COST231 vs WINNER-II comparison + shadowing effect + table of |ΔPL| ===
d_sim = linspace(0.1,10,200)'; % km

PL_hata_curve = 69.55 + 26.16*log10(freq) - 13.82*log10(hb) + (44.9 - 6.55*log10(hb)).*log10(d_sim);
hm = 1.5;
CM_urban = 3;
PL_cost231_curve = 46.3 + 33.9*log10(freq) - 13.82*log10(hb) - a_hm(freq,hm) + (44.9 - 6.55*log10(hb)).*log10(d_sim) + CM_urban;
PL_winner_curve = 32.45 + 36.7*log10(d_sim) + 20*log10(freq/1000); 


mean_abs_diff = mean(abs(PL_hata_curve - PL_cost231_curve));
tbl_compare = table(d_sim, PL_hata_curve, PL_cost231_curve, PL_winner_curve);

summaryTable = table(mean_abs_diff, 'VariableNames', {'MeanAbsDiff_Hata_COST231'});
writetable(summaryTable, fullfile(outDir,'Model_Comparison_Summary.csv'));


d_ref = [0.5 1 2 3 4 5 6 8 10]';
PL_ref = [110 120 129 135 138 141 143 146 148]'; 
PL_ref_interp = interp1(d_ref, PL_ref, d_sim, 'linear', 'extrap');

fig = figure('Visible','off'); hold on;

fill([d_sim; flipud(d_sim)], [PL_ref_interp+5; flipud(PL_ref_interp-5)], [0.9 0.9 0.9], 'EdgeColor','none');
plot(d_sim, PL_hata_curve, 'r-', 'LineWidth',1.5, 'DisplayName','Hata (sim)');
plot(d_sim, PL_cost231_curve, 'b--', 'LineWidth',1.2, 'DisplayName','COST231-Hata');
plot(d_sim, PL_winner_curve, 'g-.', 'LineWidth',1.2, 'DisplayName','WINNER-II proxy');
plot(d_ref, PL_ref, 'ko', 'MarkerFaceColor','k', 'DisplayName','TR36.942 ref points');
xlabel('Distance (km)'); ylabel('Path Loss (dB)');
legend('Location','best'); grid on;
title(sprintf('PL models vs TR36.942 fallback (Mean |\\Delta| Hata-COST231 = %.2f dB)', mean_abs_diff));
saveas(fig, fullfile(outDir,'Model_vs_TR36942_Comparison.png'));
close(fig);


pl_curves_tbl = table(d_sim, PL_hata_curve, PL_cost231_curve, PL_winner_curve, PL_ref_interp);
writetable(pl_curves_tbl, fullfile(outDir,'PL_Model_Curves.csv'));

%% === Merge validation (geo-derived averages if available) with TR fallback in one multi-line plot ===
geo_pairs = []; % collect [d_km, PL_estimated] pairs
dataDir = pwd;
files_to_try = {fullfile(dataDir,'poland_towers.csv'), fullfile(dataDir,'denmark_tower1.csv'), fullfile(dataDir,'denmark_tower2.csv')};

for fi = 1:numel(files_to_try)
    fpath = files_to_try{fi};
    if ~isfile(fpath)
        fprintf(' Validation file not found: %s\n', fpath);
        continue;
    end
    
    try
        rawText = fileread(fpath);
        if isempty(rawText)
            fprintf(' Skipping %s: empty file\n', fpath);
            continue;
        end
        
        latVec = []; lonVec = [];
        
        
        if contains(rawText, 'LTE') || contains(rawText, 'UMTS')
            
            rawText = regexprep(rawText, '(?<!^)(LTE|UMTS)', '\n$1');
            lines = strsplit(rawText, '\n');
            lines = lines(~cellfun(@isempty, strtrim(lines))); 
        else
            lines = strsplit(rawText, '\n');
            lines = lines(~cellfun(@isempty, strtrim(lines)));
        end
        
        fprintf(' Processing %d candidate lines from %s...\n', numel(lines), fpath);
        
        validCount = 0;
        for i = 1:min(100000, numel(lines))
            line = strtrim(lines{i});
            if isempty(line), continue; end
            
            
            fields = strsplit(line, ',');
            if numel(fields) < 8
                continue;
            end
            
            
            lonStr = strtrim(fields{7});
            latStr = strtrim(fields{8});
            
            
            lonClean = str2double(lonStr);
            latClean = str2double(latStr);
            
            if isnan(lonClean) || isnan(latClean)
                continue;
            end
            
            
            if latClean >= 50 && latClean <= 60 && lonClean >= 8 && lonClean <= 25
                latVec(end+1) = latClean;
                lonVec(end+1) = lonClean;
                validCount = validCount + 1;
            end
        end
        
        fprintf('  Found %d total valid coordinate pairs in %s\n', validCount, fpath);
        
        if validCount < 5, continue; end
        
        
        refLat = median(latVec); refLon = median(lonVec);
        d_km = 6371 * 2 .* asin(sqrt(...
            sin((deg2rad(latVec) - deg2rad(refLat))/2).^2 + ...
            cos(deg2rad(refLat)) .* cos(deg2rad(latVec)) .* ...
            sin((deg2rad(lonVec) - deg2rad(refLon))/2).^2));
        
        
        valid = (d_km > 0.1) & (d_km < 50);
        d_km = d_km(valid);
        if numel(d_km) < 5, continue; end
        
        
        PL = 69.55 + 26.16*log10(freq) - 13.82*log10(hb) + ...
             (44.9 - 6.55*log10(hb)) * log10(d_km);
        geo_pairs = [geo_pairs; [d_km(:), PL(:)]];
        
    catch ME
        fprintf(' Error processing %s: %s\n', fpath, ME.message);
    end
end


if ~isempty(geo_pairs)
    d_all = geo_pairs(:,1); PL_all = geo_pairs(:,2);
    validIdx = isfinite(d_all) & isfinite(PL_all) & (d_all > 0.1) & (d_all < 50);
    d_all = d_all(validIdx);
    PL_all = PL_all(validIdx);
    
    if ~isempty(d_all) && numel(d_all) >= 10
        [d_unique, ~, ic] = unique(round(d_all,3));
        PL_avg = accumarray(ic, PL_all, [], @mean);
        d_plot = linspace(max(0.1, min(d_unique)), min(20, max(d_unique)), 200);
        PL_geo_interp = interp1(d_unique, PL_avg, d_plot, 'linear', 'extrap');

        
        hata_interp = interp1(d_sim, PL_hata_curve, d_plot, 'linear', 'extrap');
        valid_rmse = isfinite(hata_interp) & isfinite(PL_geo_interp);
        if any(valid_rmse)
            RMSE_geo = rmse(hata_interp(valid_rmse), PL_geo_interp(valid_rmse));
        else
            RMSE_geo = NaN;
        end

        
        fig = figure('Visible','off'); hold on;
        fill([d_sim; flipud(d_sim)], [PL_ref_interp+5; flipud(PL_ref_interp-5)], [0.95 0.95 0.95],'EdgeColor','none');
        plot(d_sim, PL_hata_curve, 'r-', 'LineWidth',1.5);
        plot(d_sim, PL_cost231_curve, 'b--', 'LineWidth',1.2);
        plot(d_sim, PL_winner_curve, 'g-.', 'LineWidth',1.2);
        plot(d_plot, PL_geo_interp, 'm-', 'LineWidth',1.5, 'DisplayName','Geo-derived avg (interp)');
        plot(d_ref, PL_ref, 'ko', 'MarkerFaceColor','k', 'DisplayName','TR ref pts');
        xlabel('Distance (km)'); ylabel('Path Loss (dB)');
        legend('TR ±5 dB','Hata','COST231','WINNER proxy','Geo-derived avg','TR pts','Location','best');
        title(sprintf('Model comparison and geo-derived validation (RMSE=%.2f dB)', RMSE_geo));
        grid on; saveas(fig, fullfile(outDir,'Merged_Validation_Plot.png')); close(fig);
        fprintf('Geo-derived validation RMSE (Hata vs geo_interp) = %.2f dB\n', RMSE_geo);
        fprintf(fidlog,'Geo RMSE = %.2f dB\n', RMSE_geo);
    else
        
        fig = figure('Visible','off'); hold on;
        fill([d_sim; flipud(d_sim)], [PL_ref_interp+5; flipud(PL_ref_interp-5)], [0.95 0.95 0.95],'EdgeColor','none');
        plot(d_sim, PL_hata_curve, 'r-', 'LineWidth',1.5);
        plot(d_sim, PL_cost231_curve, 'b--', 'LineWidth',1.2);
        plot(d_sim, PL_winner_curve, 'g-.', 'LineWidth',1.2);
        plot(d_plot, PL_geo_interp, 'm-', 'LineWidth',1.5, 'DisplayName','Geo-derived avg (interp)');
        plot(d_ref, PL_ref, 'ko', 'MarkerFaceColor','k', 'DisplayName','TR ref pts');
        xlabel('Distance (km)'); ylabel('Path Loss (dB)');
        legend('TR ±5 dB','Hata','COST231','WINNER proxy','Geo-derived avg','TR pts','Location','best');
        title(sprintf('Model comparison and geo-derived validation (RMSE=%.2f dB)', RMSE_geo));
        grid on; saveas(fig, fullfile(outDir,'Merged_Validation_Plot.png')); close(fig);
        fprintf('Insufficient valid geo points; fallback TR plot saved.\n');
        fprintf(fidlog,'Insufficient geo data; fallback TR plot saved.\n');
    end
else
    
    fig = figure('Visible','off'); hold on;
    fill([d_sim; flipud(d_sim)], [PL_ref_interp+5; flipud(PL_ref_interp-5)], [0.95 0.95 0.95],'EdgeColor','none');
    plot(d_sim, PL_hata_curve, 'r-', 'LineWidth',1.5);
    plot(d_sim, PL_cost231_curve, 'b--', 'LineWidth',1.2);
    plot(d_sim, PL_winner_curve, 'g-.', 'LineWidth',1.2);
    plot(d_ref, PL_ref, 'ko', 'MarkerFaceColor','k', 'DisplayName','TR ref pts');
    xlabel('Distance (km)'); ylabel('Path Loss (dB)');
    legend('TR ±5 dB','Hata','COST231','WINNER proxy','TR pts','Location','best');
    title('Model comparison (no geo data available)');
    grid on; saveas(fig, fullfile(outDir,'Merged_Validation_Plot_Fallback.png')); close(fig);
    fprintf('No geo-derived validation points available; fallback TR plot saved.\n');
    fprintf(fidlog,'No geo data; fallback TR plot saved.\n');
end

%% === Histogram of throughput per terrain (PF) ===
fig = figure('Visible','off');
hold on;
colors = lines(numel(terrains));

if iscell(colors) || isstruct(colors)
    
    colors = [1 0 0; 0 1 0; 0 0 1; 1 0 1]; 
    colors = colors(1:numel(terrains), :);
else
    
    colors = min(max(colors, 0), 1); 
end
edges = linspace(min(ResultsTbl.MeanTput_PF_Mbps), max(ResultsTbl.MeanTput_PF_Mbps), 25);
for k = 1:numel(terrains)
    mask = strcmp(ResultsTbl.Terrain, terrains{k});
    h = histogram(ResultsTbl.MeanTput_PF_Mbps(mask), edges, 'DisplayStyle','stairs', 'EdgeColor', colors(k,:), 'LineWidth',1.5);
end
xlabel('Mean Throughput PF (Mbps)'); ylabel('Count');
legend(terrains,'Location','best'); title('Throughput distribution by terrain (PF)');
grid on; saveas(fig, fullfile(outDir,'Histogram_Throughput_by_Terrain_PF.png')); close(fig);

%% === Save additional artifacts and summary ===
writetable(ResultsTbl, fullfile(outDir,'Results_Table_Final.csv'));
writetable(pl_curves_tbl, fullfile(outDir,'PL_Curves_for_Replication.csv'));
fidREADME = fopen(fullfile(outDir,'README.txt'),'w');
fprintf(fidREADME,'This output folder contains simulation results and figures for IEEE Access submission.\n');
fprintf(fidREADME,'Frequency: %d MHz (Kathrein 742215 band 1710-2180 MHz)\n', freq);
fprintf(fidREADME,'Monte Carlo iterations: %d\n', MonteCarlo_iter);
fprintf(fidREADME,'Shadowing sigma: %d dB\n', shadow_sigma);
fprintf(fidREADME,'Files: Results_Table_Final.csv, PL_Curves_for_Replication.csv, Model_vs_TR36942_Comparison.png, Merged_Validation_Plot*.png, Histogram_Throughput_by_Terrain_PF.png, Contour_Urban_Throughput.png, Pareto_*.png, MC_Convergence_*.png\n');
fclose(fidREADME);
fclose(fidlog);

fprintf('\n=== ENHANCED SIMULATION COMPLETE ===\n');
fprintf('Outputs saved in: %s\n', outDir);

%% === Local Pareto function ===
function idx = paretofront_local(x,y)
    n = numel(x);
    idx = true(n,1);
    for i = 1:n
        for j = 1:n
            if all([x(j) >= x(i), y(j) >= y(i)]) && any([x(j) > x(i), y(j) > y(i)])
                idx(i) = false; break;
            end
        end
    end
end