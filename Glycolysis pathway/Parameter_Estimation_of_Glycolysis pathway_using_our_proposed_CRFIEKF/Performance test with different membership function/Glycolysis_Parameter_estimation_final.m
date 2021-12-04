%%%
% Script version of the Regularized Fuzzy Inferred Extended Kalman Filter for parameter estimation of glycolysis pathway
%   - Number of States: 5
%   - Number of estimated prameters: 47
%   - Number of input to FIS: 9
%   - Number of output from FIS: 6
%
% works!!

%%% Initialization

clear all;
close all;
clc;

% model parameters
m2	  =	0.7152	;	% glc_6P
m3	  =	0.6143	;	% glc arbitary
m4	  =	0.922	;	% f_6P
m5	  =	0.6691	;	% f_16_BP 
m6	  =	0.1878	;	% f_26_BP arbitary
m7	  =	0.3506	;	% dhap arbitary
m8	  =	0.5922	;	% ga_3P
m9	  =	0.9618	;	% BPG_13
m10	  =	0.2324;%0.9684	;	% PG_3
m11	  =	0.2419	;	% PG_2
m12	  =	0.24725;%0.9735	;	% PEP
m13	  =	0.9615	;	% PYR arbitary
m14	  =	0.5368	;	% LACTATE arbitary (only produced)

m31	  =	0.7354	;	% ATP
m32	  =	0.5286	;	% ADP arbitary
m33	  =	0.3492	;	% NADH
m34	  =	0.1416	;	% NAD


%%%%% ENZYME/GENE INITIAL CONCENTRATION g in nM (Assumed all) %%%%%
g3	=	0.735;	%Hexokinase				
g4	=	0.631;	%Glucose-6-phosphatase				
g5	=	0.34;	%Phosphoglucoisomerase				
g6	=	0.6275;	%Phosphofructokinase_1				
g7	=	0.784;	%Fructose-1,6-bisphosphatase				
g8	=	0.322;	%Phosphofructokinase_2				
g9	=	0.503;	%Fructose-2,6-bisphosphatase			
g10	=	0.719;	%Aldolase
g11	=	0.281;	%Triose_phosphate_isomerase					
g12	=	0.7627;	%Glyceraldehyde-3-phosphate_dehydrogenase				
g13	=	0.4552;	%Phosphoglycerate_kinase				
g14	=	0.105;	%Phosphoglycerate_mutase				
g15	=	0.125;	%Enolase					
g16	=	0.6981;	%Pyruvate_kinase				
g17	=	0.1176;	%Lactate_dehydrogenase	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% METABOLIC REACTION SUBSTRATE CONSTANT km in nM %%%%%

km3	 =	0.052076	;
km4	 =	0.4479	;
km5	 =	0.052373	;
km6	 =	0.4289	;
km7	 =	0.054248	;
km8	 =	0.4654	;
km9	 =	0.053261	;
km10	 =	0.07412	;
km11	 =	0.1	;
km12	 =	0.4415	;
km13	 =	0.09692	;
km14	 =	0.7	;
km15	 =	0.051036	;
km16	 =	0.4729	;
km17	 =	0.2;
km18	 =	0.4658	;
km19	 =	0.053377	;
km20	 =	0.3951	;
km21	 =	0.073413	;
km22	 =	0.42947	;
km23	 =	0.13745	;
km24	 =	0.7	;
 
% 
% %%%%% METABOLIC REACTION ENZYME RATE CONATANT K in S^-1  %%%%%
% 
 
K3	 =	0.02	;	% hk
K4	 =	0.05	;	% g6Pase
K5	 =	0.0812	;	% pgi
K6	 =	0.09	;	% pfk1
K7	 =	0;%0.017	;	% f16Bpase
K8	 =	0.09	;	% pfk2
K9	 =	0;%0.0766	;	% f26Bpase
K10	 =	0.0826	;	% ald
K11	 =	0.0899	;	% tpi
K12	 =	0.01	;	% gcld3PD
K13	 =	0.0532	;	% pglc_kn
K14	 =	0.06	;	% pglc_m
K15	 =	0.0727	;	% enl
K16	 =	0.2	;	% pyrk
K17	 =	0.01	;	% lacd


F4	=	0.5;
F5	=	0.5;
F7	=	0.5;
F8	=	0.5;
F11	=	0.5;
F15	=	0.5;
F16	=	0.5;
F19	=	0.5;
F20	=	0.5;
F23	=	0.5;
F27	=	0.5;
F28	=	0.5;


%%% Number of parameters = 73 
p = [m3,m4,m6,m7,m8,m10,m11,m12,m14,m31,m34,g3,g4,g5,g6,g7,g8,g9,g10,g11,g12,g13,g14,g15,g16,g17,km3,km4,km5,km6,km7,km8,km9,km10,km11,km12,km13,km14,km15,km16,km17,km18,km19,km20,km21,km22,km23,km24,K3,K4,K5,K6,K8,K10,K11,K12,K13,K14,K15,K16,K17,F4,F5,F7,F8,F11,F15,F16,F19,F20,F23,F27,F28];

% Model used for estimation
model =@(t,y) Glycolysis_extended_model_final(t,y);
jac=@(t,y) Glycolysis_model_jacobean_final(t,y);


% Set up time vector 
T  = 10;
Ts = 0.002;
Tspan = [0 Ts];
time = [0:Ts:T]';

% Find the initial conditions (equilibrium point)
x0 = [0.7152, 0.6143, 0.922, 0.6691, 0.1878, 0.3506];


% Dimensions
S=length(x0);
P=length(p);
N=S+P;


%%% Obtain measurement signal

% Output matrix

H = zeros(6,79);

H(1,12) = 1;
H(2,15) = 1;
H(3,1) = 1;
H(4,8) = 1;
H(5,16) = 1;
H(6,4) = 1;

   
% Measurement noise covariance
R = eye(6)*0.001;
        
% Process noise
Q=0.01*eye(79);

% Initial conditions
xe0 = [0.7152	0.6143	0.922	0.6691	0.1878	0.3506	0.5922	0.9618	0.2324	0.2419	0.24725	0.9615	0.5368	0.7354	0.5286	0.3492	0.1416	0.735	0.631	0.34	0.6275	0.784	0.322	0.503	0.719	0.281	0.7627	0.4552	0.105	0.125	0.6981	0.1176	0.09	"	0.1"	0.1	0.9	0.1	0.15	0.5	0.01	0.01	0.8	0.01	0.9	0.1	0.4	0.01	0.8	0.3	0.01	0.5	0.1	0.9	0.1	0.5	0.7	0.5	0.5	0.5	0.1	0.4	0.3	0.1	0.2	0.2	0.1	0.09	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5]';
Pe0=eye(79);


% State constraints Dx<=d
D = -eye(N);
d = -0.1*ones(N,1);

alpha = 0.56;%rand(1);
reg_mat = Pe0;

disp('Parameter estimation started...please wait...')
[xe_Gaussian, estimated_parameter_Gaussian, Error]= CRFIEKF_with_Gaussian_FIS(time, model, jac, H , Q, R, xe0, Pe0, D, d, alpha, reg_mat);
disp('Parameter estimation complited with Gaussian membership function')
[xe_Gen_Bell, estimated_parameter_Gen_Bell, Error]= CRFIEKF_with_Gen_Bell_FIS(time, model, jac, H , Q, R, xe0, Pe0, D, d, alpha, reg_mat);
disp('Parameter estimation complited with Generalized Bell membership function')
[xe_Triangular, estimated_parameter_Triangular, Error]= CRFIEKF_with_Triangular_FIS(time, model, jac, H , Q, R, xe0, Pe0, D, d, alpha, reg_mat);
disp('Parameter estimation complited with Traiangular membership function')
[xe_Trapiziodal, estimated_parameter_Trapiziodal, Error]= CRFIEKF_with_Trapiziodal_FIS(time, model, jac, H , Q, R, xe0, Pe0, D, d, alpha, reg_mat);
disp('Parameter estimation complited with Trapiziodal membership function')

disp('Parameter estimation completed');
disp('Validation with Hypoxia WET lab data started... please waight...');

% load('Hypoxia_WET_lab_validation_02_10_2018.mat')
x_v_fis_Gaussian = Hypoxia_simlation(estimated_parameter_Gaussian);
x_v_fis_Gen_Bell = Hypoxia_simlation(estimated_parameter_Gen_Bell);
x_v_fis_Triangular = Hypoxia_simlation(estimated_parameter_Triangular);
x_v_Trapiziodal = Hypoxia_simlation(estimated_parameter_Trapiziodal);
disp('Validation with Hypoxia WET lab data completed');
disp('Generating plot...');

Parameter_estimation_FIS_membership_vs_Kinoshita_plot(x_v_fis_Gaussian, x_v_fis_Gen_Bell, x_v_fis_Triangular, x_v_Trapiziodal);
% validation_plot();



