function [xe, M, error] = CRFIEKF_Tikhonov_FIS_4_output(time, model, jac, Hk , Q, R, xe0, Pe0, D, d, alpha)

% Detect model size
N=size(Q,2);

% Set initial value
xe=zeros(N, 1);


xe(:,1)=xe0;
constrained_xe = xe;

Pe=cell(1,5000);
I=eye(N);
Pe{1}=Pe0;


%%% Create Fuzzy inference system
fis_glc_pyrk = readfis('fuzzy_rule_glc_pyrk');
options = evalfisOptions;
options.NumSamplePoints = 200;

% fis_hk_atp_glc = readfis('fuzzy_rule_hk_atp_glc');
% options = evalfisOptions;
% options.NumSamplePoints = 200;


fis_gap_g3pdh_nad = readfis('fuzzy_rule_gap_g3pdh_nad');
options = evalfisOptions;
options.NumSamplePoints = 200;


fis_f6p_pfk1 = readfis('fuzzy_rule_f6p_pfk1');
options = evalfisOptions;
options.NumSamplePoints = 200;


observed = zeros(4,1);
% u(1) = 0;
% u(1) = 0.5; 
% u(2) = 0.5; 
% u(3) = 0.01; %((k30*x16)/(k31*x3) - 1)/F4; 
% u(4) = 0.01; %((k45*x23)-(k47*x7))/k46; 
error  = 0;

for ii=2:5000
    % **1** Compute linearization around previous a posteriori estimate
    Fk=jac(time(ii-1),xe(:,ii-1));
    
    % **2** A priori estimate of the current state
    opts=odeset('RelTol',1e-8,'AbsTol',1e-10);
    sol=ode113(model,[time(ii-1) time(ii)], xe(:,ii-1),opts);
    S=size(sol.y);
    xp=sol.y(:,S(2));
    
    % **3** Restrict a priori estimate within range [0,1]
    for j = 1:79
        if (xp(j, 1))<0.1
            xp(j, 1) = 0.1;
        end
    end
   
    for j = 1:79
        if (xp(j, 1))>1.0
            xp(j, 1) = 1.0;
        end   
    end

    
    % **4** Integrates the continuos-time differential Lyapunov equation
    Pp=integration_cov_mat(Fk,Q,[time(ii-1) time(ii)],Pe{ii-1});
    Pp=Pp{numel(Pp)};
    Pp=(Pp+Pp')/2; % Symmetrize
    
    % **5** Compute optimal Kalman gain
    Lk=Pp*Hk'*(Hk*Pp*Hk'+R)^(-1);
    
    
    % **6** Generate measument from Fuzzy inference system
% %     observed(1,ii-1) = evalfis([xp(2,1) xp(23,1)],fis_glc_pyrk,options);
% %     actual = evalfis([xp(18,1) xp(14,1) xp(2,1) ],fis_hk_atp_glc,options);
% %     observed(2:3,ii-1) = [actual(2); actual(1)]; 
% %     observed(4:5,ii-1) = evalfis([xp(7,1) xp(21,1) xp(17,1)],fis_gap_g3pdh_nad,options);
% %     observed(6,ii-1) = evalfis([xp(3,1) xp(19,1) xp(14,1)],fis_f6p_pfk1,options);


    observed(1,ii-1) = evalfis([xp(2,1) xp(23,1)],fis_glc_pyrk,options);
    observed(2:3,ii-1) = evalfis([xp(7,1) xp(21,1) xp(17,1)],fis_gap_g3pdh_nad,options);
    observed(4,ii-1) = evalfis([xp(3,1) xp(19,1) xp(14,1)],fis_f6p_pfk1,options);

    
    
    % **7** Incorporate fuzzy measurement to get a posteriori estimate of current prediction 
    xe(:,ii)=xp+Lk*(observed(1:4,ii-1) - Hk*xp);
    Pe{ii}=(I-Lk*Hk)*Pp*(I-Lk*Hk)'+(Lk*R*Lk');
    Pe{ii}=(Pe{ii}+Pe{ii}')/2; % Symmetrize

    
    % **8** If condition number > Threshold ==> replace xe with Regularized xe
    reg_mat = pinv(Pe{ii});
    N_k = Hk'*(R^-1)*Hk + Pe{ii};
    n = norm(N_k);
    con(ii) = norm(N_k) * norm(pinv(N_k));
    
    if (con(ii) > 9.30e3)
        test =  (Hk'*R^(-1)*Hk + Pe{ii} + alpha*reg_mat)^(-1);
        K_trkf = test  * (Hk'*R^(-1));
  
        x_trkf = xp + K_trkf*(observed(1:4,ii-1) - Hk*xp);
  
        P_trkf = (I - K_trkf*Hk)*Pe{ii};
  
        xe(:,ii) = x_trkf(:,1);
        Pe{ii} = P_trkf;
        Pe{ii}=(Pe{ii}+Pe{ii}')/2; 
    end
    
    
    % **9** If constraints not satisfied ==> replace Regularized xe with constrained regularized xe   
    if ~all(D*xe(:,ii)<=d)
        % Check invertibility and positive definiteness of Pe{ii}
        if all(eig(Pe{ii})>zeros(size(Pe{ii},1),1))
            W=inv(Pe{ii});
            W=(W+W')/2;
        else
            disp('Warning. Non positive-definite Pe.')
            W=pinv(Pe{ii});
            W=(W+W')/2;   
            break;
        end
        % Solve the quadratic programming problem
        oo = optimset('Algorithm','interior-point-convex','Display','off');
        z = quadprog(2*W, zeros(size(xe(:,ii-1))), D, d-D*xe(:,ii),[],[],[],[],[],oo);
        % Replace unconstrained estimate with constrained estimate

        constrained_xe(:,ii) = xe(:,ii) + z;
        xe(:,ii) = constrained_xe(:,ii);
    else
        xe(:,ii) = xe(:,ii);
    end  
    
    % **10** Restrict upper bound of constrained regularized xe to 1.
    for j = 1:79
        if (xe(j, ii))>1.0
            xe(j, ii) = 1.0;
        end   
    end
       

    % **11** Compute RMS error
    
% %     actual(1:6,ii-1) = [xe(12,ii); xe(15,ii); xe(1,ii); xe(8,ii); xe(16,ii); xe(4,ii)];
% %     
% %     sse = 0;
% %     for k = 1:6
% %        sse = sse +  (observed(k,ii-1) - actual(k,ii-1))^2;
% %     end
% %     
% %     mse = sse/6;
% %     rmse = sqrt(mse);    
% %     error = error + rmse;    


    actual(1:4,ii-1) = [xe(12,ii); xe(8,ii); xe(16,ii); xe(4,ii)];
    
    sse = 0;
    for k = 1:4
       sse = sse +  (observed(k,ii-1) - actual(k,ii-1))^2;
    end
    
    mse = sse/4;
    rmse = sqrt(mse);    
    error = error + rmse;  
    
    if(mod(ii, 1000) == 0)
        fprintf("Iteration %d completed\n",ii);
    end
end

% fprintf("\n\n parameter estimation process ends...\n\n");
%%% Estimation results

for i= 1:47
   M(i,1) = mean(xe(i+32,ii-200:ii-1)); 
end

end