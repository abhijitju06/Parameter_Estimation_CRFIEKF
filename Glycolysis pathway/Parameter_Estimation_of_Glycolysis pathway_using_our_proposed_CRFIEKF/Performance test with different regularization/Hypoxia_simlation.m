function [x_v] = Hypoxia_simlation(k)

% clc
% clear
% close all

warning('off','all');

global km3	km4	km5	km6	km7	km8	km9	km10	km11	km12	km13	km14	km15	km16	km17	km18	km19	km20	km21	km22	km23	km24	K3	K4	K5	K6	K7	K8	K9	K10	K11	K12	K13	K14	K15	K16	K17	F4	F5	F7	F8	F11	F15	F16	F19	F20	F23	F27	F28;

km3=	k(1);
km4=	k(2);
km5=	k(3);
km6=	k(4);
km7=	k(5);
km8=	k(6);
km9=	k(7);
km10=	k(8);
km11=	k(9);
km12=	k(10);
km13=	k(11);
km14=	k(12);
km15=	k(13);
km16=	k(14);
km17=	k(15);
km18=	k(16);
km19=	k(17);
km20=	k(18);
km21=	k(19);
km22=	k(20);
km23=	k(21);
km24=	k(22);
K3=     k(23);
K4=     k(24);
K5=     k(25);
K6=     k(26);
K7=     0;
K8=     k(27);
K9=     0;
K10=	k(28);
K11=	k(29);
K12=	k(30);
K13=	k(31);
K14=	k(32);
K15=	k(33);
K16=	k(34);
K17=	k(35);
F4=     k(36);
F5=     k(37);
F7=     k(38);
F8=     k(39);
F11=	k(40);
F15=	k(41);
F16=	k(42);
F19=	k(43);
F20=	k(44);
F23=	k(45);
F27=	k(46);
F28=	k(47);










%%%%% METABOLITES INITIAL CONCENTRATION in nM %%%%%

m2	  =	0.6152	;	% glc_6P
m3	  =	0.6143	;	% glc arbitary
m4	  =	0.822	;	% f_6P
m5	  =	0.4691	;	% f_16_BP 
m6	  =	0.1878	;	% f_26_BP arbitary
m7	  =	0.4506	;	% dhap arbitary
m8	  =	0.5922	;	% ga_3P
m9	  =	0.9618	;	% BPG_13
m10	  =	0.29324;    %0.9684	;	% PG_3
m11	  =	0.2419	;	% PG_2
m12	  =	0.64725;    %0.9735	;	% PEP
m13	  =	0.657615;	% PYR arbitary
m14	  =	0.5368	;	% LACTATE arbitary (only produced)

m31	  =	0.5354	;	% ATP
m32	  =	0.3286	;	% ADP arbitary
m33	  =	0.3492	;	% NADH
m34	  =	0.1416	;	% NAD



%%%%% ENZYME/GENE INITIAL CONCENTRATION g in nM (Assumed all) %%%%%


g3	 =	0.735;      %Hexokinase
g4	=	0.4631;     %Glucose-6-phosphatase				
g5	=	0.2534;     %Phosphoglucoisomerase				
g6	 =	0.74;       %Phosphofructokinase_1
g7	=	0.784;      %Fructose-1,6-bisphosphatase				
g8	=	0.322;      %Phosphofructokinase_2				
g9	=	0.503;      %Fructose-2,6-bisphosphatase	
g10	 =	0.719;      %Aldolase arbitrary 0.809 actual
g11	=	0.281;      %Triose_phosphate_isomerase				
g12	 =	0.7627;     %Glyceraldehyde-3-phosphate_dehydrogenase
g13	 =	0.64552;	%Phosphoglycerate_kinase
g14	=	0.105;      %Phosphoglycerate_mutase				
g15	=	0.125;      %Enolase				
g16	 =	0.6981;     %Pyruvate_kinase
g17	 =	0.1176;     %Lactate_dehydrogenase

x_v = zeros(1,32);
x_v(1:17) = [m2 m3 m4 m5 m6 m7 m8 m9 m10 m11 m12 m13 m14 m31 m32 m33 m34];
x_v(18:32) = [g3 g4 g5 g6 g7 g8 g9 g10 g11 g12 g13 g14 g15 g16 g17];

% uu = zeros(1, 27); 

T  = 10;
Ts = 0.002;
Tspan = [0 Ts];
t = [0:Ts:T]';
data = xlsread('Activities_hypoxia_normalized.xlsx');
tic

for ii = 1:(T/Ts)

x_v(ii,18) = data(ii,1);
x_v(ii,21) = data(ii,2);
x_v(ii,25) = data(ii,3);
x_v(ii,27) = data(ii,4);
x_v(ii,28) = data(ii,5);
x_v(ii,31) = data(ii,6);
x_v(ii,32) = data(ii,7); 

                                                            
[tsim, xsim] = ode23tb('Glycolysis_hypoxia_model', Tspan, x_v(ii,:));
x_v(ii+1, :) = real(xsim(end, :));

for j = 1:32
   if (x_v(ii+1, j))<0.1
      x_v(ii+1, j) = 0.1;
   end
end
   
for j = 1:32
   if (x_v(ii+1, j))>1.0
      x_v(ii+1, j) = 1.0;
   end   
end
%  if (ii<(T/Ts))
%         fprintf('Iteration %d is completed...\n\n',ii);
%     else
%         fprintf('Iteration %d is completed and finally simulation is completed.\n\n',ii);
%  end

end

% fprintf('Hypoxia simulation completed...\n');


% load('Hypoxia_WET_lab_validation_02_10_2018.mat')
% 
% validation_plot();

end
