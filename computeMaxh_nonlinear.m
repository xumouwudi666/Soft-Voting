function [maxh_matrix, AngleALL] = computeMaxh_nonlinear(clstack, corrstack, K, p, method, sigma, a)
% Nonlinear weighted method for computing angle between equivalent lines
% method: 'gaussian' or 'sigmoid', specifies the weighting method
% sigma: for gaussian weighting, controls the width of the gaussian function
% threshold: the threshold for correlation to apply non-zero weights
% a: parameter for sigmoid function

% Initialize matrices
maxh_matrix = zeros(K, K);
AngleALL = zeros(K, K);
weight = ones(K, K, K, 'double'); % Initialize weights to 1
T = 60;
tho = 180 / T;
count=0;

% Compute initial maxh_matrix using original method
maxh_matrix_old = computeMaxhandSts_ang(clstack, corrstack, K);
%load('maxh_matrix_old_SNR8.mat');

% ========== 新增：迭代控制参数 ==========
max_iter = 30;                    % 最大迭代次数
patience = 5;                     % 振荡检测窗口
relative_std_threshold = 0.01;    % 相对标准差阈值（1%）
error_history = [];               % 记录误差历史
% =======================================
%原先的
Thresh = prctile(maxh_matrix_old, p, 'all');



flagW = 1;
sum_countW_old=0;
sum_maxh_matrix_old = sum(sum(maxh_matrix_old));
% Iteration for updating weights and computing angles
while count < max_iter
    sum_countW = 0;
    
    for k1 = 1:K-1
        for k2 = k1+1:K
            h = zeros(1, T);
            angles = -1 * ones(1, K);
            countW = 0;
            

            
            
            for k3 = 1:K
                if k3 == k1 || k3 == k2
                    continue;
                end
                
                % Update weight based on the chosen method
                %if method == 'gaussian'
                     %weight(k1, k2, k3) = exp(- (corrstack(k1, k3) - Thresh)^2 / (2 * sigma^2));
                % elseif method == 'sigmoid'
                     weight(k1, k2, k3) = 1 / (1 + exp(-a * (corrstack(k1, k3) - Thresh)));
                     %weight(k1, k2, k3) = 1/(1+exp(-a*(maxh_matrix_old(k1,k3) - Thresh)));
                % end
                %if method == 'tanh'
                      %weight(k1, k2, k3) = 0.5 * (1 + tanh(a * (corrstack(k1, k3) - Thresh)));
                      %weight(k1, k2, k3) = 0.5 * (1 + tanh(a * (corrstack(k1, k3) - Thresh_corr)));
                % elseif strcmp(method, 'power')
                 %if corrstack(k1, k3) >= Thresh
                    %normalized_corr = (corrstack(k1, k3) - Thresh) / (1 - Thresh);
                    %weight(k1, k2, k3) = normalized_corr^a;  % a在power方法中代表幂次
                 %else
                    %weight(k1, k2, k3) = 0;
                 %end
                
                if maxh_matrix_old(k1, k3) > Thresh && maxh_matrix_old(k2, k3) > Thresh
                    countW = countW + 1;
                    
             
                    
                    
                else
                    weight(k1, k2, k3) = 0;
                end
                     %countW = countW + 1;
                
                
                % Compute angle using current weight
                [angle12, flag] = ComputeAngleCorr(clstack(k1, k2), clstack(k2, k1), corrstack(k1, k2), clstack(k1, k3), clstack(k3, k1), corrstack(k1, k3), clstack(k2, k3), clstack(k3, k2), corrstack(k2, k3));
                angles(k3) = angle12;
                
                if flag
                    ang = (1:T) * 180.0 / T;
                    h = h + double(weight(k1, k2, k3)) * exp(-(ang - angle12).^2 / (2 * tho * tho)); % Nonlinear weighted voting
                end
            end
            
            %if(countW == 0)
                %flagW = 0;
                %break;
            if(countW == 0)
                maxh_matrix(k1, k2) = 0;  
                maxh_matrix(k2, k1) = 0;
                continue;
            else
                
                h = h./countW; % Modified by Lu
                sum_countW = sum_countW+countW;
                
            end

            [maxh,idx] = max(h);
            
            
            ang = idx*tho - tho/2;
            maxh_matrix(k1, k2) = maxh;
            maxh_matrix(k2, k1) = maxh;
            
            AngleALL(k1, k2) = ang;
            AngleALL(k2, k1) = ang;

        end
        
        %if(flagW == 0)
           % break;
        %end

    end
    
      % ========== 新增：振荡检测 ==========
    current_error = max(max(abs(maxh_matrix - maxh_matrix_old)));
    error_history = [error_history; current_error];
    
    % 检查振荡（只在积累够数据后）
    oscillation_detected = false;
    if length(error_history) >= patience
        recent_errors = error_history(end-patience+1:end);
        recent_mean = mean(recent_errors);
        recent_std = std(recent_errors);
        relative_std = recent_std / recent_mean;
        
        if relative_std < relative_std_threshold
            fprintf('振荡检测: 最近%d次误差mean=%.6f, std=%.6f, 相对std=%.2f%% < %.2f%%\n', ...
                    patience, recent_mean, recent_std, relative_std*100, relative_std_threshold*100);
            oscillation_detected = true;
        end
    end
    % =====================================
    
    %fprintf('max(maxh_matrix_old)==%f,   max(maxh_matrix)== %f,   error==%f\n',max(max(maxh_matrix_old)), max(max(maxh_matrix)), max(max(abs(maxh_matrix - maxh_matrix_old))));
    % =====================================
    
    fprintf('max(maxh_matrix_old)==%f,   max(maxh_matrix)== %f,   error==%f\n', ...
            max(max(maxh_matrix_old)), max(max(maxh_matrix)), current_error);
    
    % 修改：添加振荡检测到终止条件
    %if  (flagW == 0) || oscillation_detected
    if  (current_error < 0.01) ||(flagW == 0) || oscillation_detected
        fprintf('终止原因: error<0.06=%d, flagW=0=%d, oscillation=%d\n', ...
            current_error<0.06, flagW==0, oscillation_detected);
        break;
    %if(max(max(abs(maxh_matrix - maxh_matrix_old)))<0.06) || (flagW == 0)
    %if (sum_countW_old == sum_countW)
       % break;
    else
        sum_maxh_matrix = sum(sum(maxh_matrix));
        fprintf('sum_countW_old==%f,   sum_countW==%f\n',sum_countW_old,sum_countW);
        fprintf('sum_maxh_matrix_old== %f,   sum_maxh_matrix==%f\n',sum_maxh_matrix_old, sum_maxh_matrix);
        maxh_matrix_old = maxh_matrix;
        sum_maxh_matrix_old = sum_maxh_matrix;
        sum_countW_old = sum_countW;
        count= count+1;
        fprintf('iter==%f\n',count);
        %figure; imagesc(maxh_matrix);
        %figure; histogram(maxh_matrix);
    end
end
% 新增：最终输出信息
fprintf('========== 迭代完成 ==========\n');
fprintf('总迭代次数: %d\n', count);
fprintf('最终误差: %.6f\n', current_error);
if  oscillation_detected
    fprintf('终止原因: 检测到稳定振荡\n');
elseif count >= max_iter
    fprintf('终止原因: 达到最大迭代次数\n');
end
fprintf('==============================\n');

end
