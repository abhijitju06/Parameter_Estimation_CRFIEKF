function jac = Glycolysis_model_jacobean(t,y)


m2	=	y(1);     %	glc_6P					
m3	=	y(2);     %	glc				
m4	=	y(3);     %	f_6P					
m5	=	y(4);     %	f_16_BP					
m6	=	y(5);     %	f_26_BP				
m7	=	y(6);     %	dhap				
m8	=	y(7);     %	ga_3P					
m9	=	y(8);     %	BPG_13					
m10	=	y(9);	%	PG_3					
m11	=	y(10);	%	PG_2					
m12	=	y(11);	%	PEP					
m13	=	y(12);	%	PYR				
m14	=	y(13);	%	LACTATE		
				
m31	=	y(14);	%	ATP					
m32	=	y(15);	%	ADP			
m33	=	y(16);	%	NADH					
m34	=	y(17);	%	NAD	



%%%%% ENZYME/GENE INITIAL CONCENTRATION g in nM  %%%%%
			
g3	=	y(18);	%Hexokinase				
g4	=	y(19);	%Glucose-6-phosphatase				
g5	=	y(20);	%Phosphoglucoisomerase				
g6	=	y(21);	%Phosphofructokinase_1				
g7	=	y(22);	%Fructose-1,6-bisphosphatase				
g8	=	y(23);	%Phosphofructokinase_2				
g9	=	y(24);	%Fructose-2,6-bisphosphatase				
g10	=	y(25);	%Aldolase
g11	=	y(26);	%Triose_phosphate_isomerase				
g12	=	y(27);	%Glyceraldehyde-3-phosphate_dehydrogenase				
g13	=	y(28);	%Phosphoglycerate_kinase				
g14	=	y(29);	%Phosphoglycerate_mutase				
g15	=	y(30);	%Enolase				
g16	=	y(31);	%Pyruvate_kinase				
g17	=	y(32);	%Lactate_dehydrogenase	

%%%%% Michaelis-Menten constant %%%%%

km3	 =	y(33);
km4	 =	y(34);
km5	 =	y(35);
km6	 =	y(36);
km7	 =	y(37);
km8	 =	y(38);
km9	 =	y(39);
km10 =	y(40);
km11 =	y(41);
km12 =	y(42);
km13 =	y(43);
km14 =	y(44);
km15 =	y(45);
km16 =	y(46);
km17 =	y(47);
km18 =	y(48);
km19 =	y(49);
km20 =	y(50);
km21 =	y(51);
km22 =	y(52);
km23 =	y(53);
km24 =	y(54);
 

%%%%% METABOLIC REACTION ENZYME RATE CONATANT K in S^-1  %%%%%


K3	 =	y(55);	% hk
K4	 =	y(56);	% g6Pase
K5	 =	y(57);	% pgi
K6	 =	y(58);	% pfk1
K7	 =	0;%0.017	;	% f16Bpase
K8	 =	y(59);	% pfk2
K9	 =	0;%0.0766	;	% f26Bpase
K10	 =	y(60);	% ald
K11	 =	y(61);	% tpi
K12	 =	y(62);	% gcld3PD
K13	 =	y(63);	% pglc_kn
K14	 =	y(64);	% pglc_m
K15	 =	y(65);	% enl
K16	 =	y(66);	% pyrk
K17	 =	y(67);	% lacd

%%%%% FEEDBACK CONATANT %%%%%

F4	=	y(68);
F5	=	y(69);
F7	=	y(70);
F8	=	y(71);
F11	=	y(72);
F15	=	y(73);
F16	=	y(74);
F19	=	y(75);
F20	=	y(76);
F23	=	y(77);
F27	=	y(78);
F28	=	y(79);


% Preallocate Jacobian
jac=zeros(79,79);
% Set nonzero elements



jac(1, 1) = (K4*g4*m2*(F7*m2 + 1))/(km4 + m2)^2 - (K5*g5)/((km5 + m2)*(F8*m5 + 1)) - (F7*K4*g4*m2)/(km4 + m2) - (K4*g4*(F7*m2 + 1))/(km4 + m2) + (K5*g5*m2)/((km5 + m2)^2*(F8*m5 + 1)) - (F5*K3*g3*m3*m31*(F4*m32 + 1))/((km3 + m3*m31)*(F5*m2 + 1)^2);

jac(1, 2) = (K3*g3*km3*m31*(F4*m32 + 1))/((km3 + m3*m31)^2*(F5*m2 + 1));

jac(1, 3) = (K5*g5*km6)/((km6 + m4)^2*(F11*m5 + 1));

jac(1, 4) = (F8*K5*g5*m2)/((km5 + m2)*(F8*m5 + 1)^2) - (F11*K5*g5*m4)/((km6 + m4)*(F11*m5 + 1)^2);

jac(1, 14) = (K3*g3*km3*m3*(F4*m32 + 1))/((km3 + m3*m31)^2*(F5*m2 + 1));

jac(1, 15) = (F4*K3*g3*m3*m31)/((km3 + m3*m31)*(F5*m2 + 1));

jac(1, 18) = (K3*m3*m31*(F4*m32 + 1))/((km3 + m3*m31)*(F5*m2 + 1));

jac(1, 19) = -(K4*m2*(F7*m2 + 1))/(km4 + m2);

jac(1, 20) = (K5*m4)/((km6 + m4)*(F11*m5 + 1)) - (K5*m2)/((km5 + m2)*(F8*m5 + 1));

jac(1, 33) = -(K3*g3*m3*m31*(F4*m32 + 1))/((km3 + m3*m31)^2*(F5*m2 + 1));

jac(1, 34) = (K4*g4*m2*(F7*m2 + 1))/(km4 + m2)^2;

jac(1, 35) = (K5*g5*m2)/((km5 + m2)^2*(F8*m5 + 1));

jac(1, 36) = -(K5*g5*m4)/((km6 + m4)^2*(F11*m5 + 1));

jac(1, 55) = (g3*m3*m31*(F4*m32 + 1))/((km3 + m3*m31)*(F5*m2 + 1));

jac(1, 56) = -(g4*m2*(F7*m2 + 1))/(km4 + m2);

jac(1, 57) = (g5*m4)/((km6 + m4)*(F11*m5 + 1)) - (g5*m2)/((km5 + m2)*(F8*m5 + 1));

jac(1, 68) = (K3*g3*m3*m31*m32)/((km3 + m3*m31)*(F5*m2 + 1));

jac(1, 69) = -(K3*g3*m2*m3*m31*(F4*m32 + 1))/((km3 + m3*m31)*(F5*m2 + 1)^2);

jac(1, 70) = -(K4*g4*m2^2)/(km4 + m2);

jac(1, 71) = (K5*g5*m2*m5)/((km5 + m2)*(F8*m5 + 1)^2);

jac(1, 72) = -(K5*g5*m4*m5)/((km6 + m4)*(F11*m5 + 1)^2);

jac(2, 1) = (K4*g4*(F7*m2 + 1))/(km4 + m2) + (F7*K4*g4*m2)/(km4 + m2) - (K4*g4*m2*(F7*m2 + 1))/(km4 + m2)^2 + (F5*K3*g3*m3*m31*(F4*m32 + 1))/((km3 + m3*m31)*(F5*m2 + 1)^2);

jac(2, 2) = -(K3*g3*km3*m31*(F4*m32 + 1))/((km3 + m3*m31)^2*(F5*m2 + 1));

jac(2, 14) = -(K3*g3*km3*m3*(F4*m32 + 1))/((km3 + m3*m31)^2*(F5*m2 + 1));

jac(2, 15) = -(F4*K3*g3*m3*m31)/((km3 + m3*m31)*(F5*m2 + 1));

jac(2, 18) = -(K3*m3*m31*(F4*m32 + 1))/((km3 + m3*m31)*(F5*m2 + 1));

jac(2, 19) = (K4*m2*(F7*m2 + 1))/(km4 + m2);

jac(2, 33) = (K3*g3*m3*m31*(F4*m32 + 1))/((km3 + m3*m31)^2*(F5*m2 + 1));

jac(2, 34) = -(K4*g4*m2*(F7*m2 + 1))/(km4 + m2)^2;

jac(2, 55) = -(g3*m3*m31*(F4*m32 + 1))/((km3 + m3*m31)*(F5*m2 + 1));

jac(2, 56) = (g4*m2*(F7*m2 + 1))/(km4 + m2);

jac(2, 68) = -(K3*g3*m3*m31*m32)/((km3 + m3*m31)*(F5*m2 + 1));

jac(2, 69) = (K3*g3*m2*m3*m31*(F4*m32 + 1))/((km3 + m3*m31)*(F5*m2 + 1)^2);

jac(2, 70) = (K4*g4*m2^2)/(km4 + m2);

jac(3, 1) = (K5*g5*km5)/((km5 + m2)^2*(F8*m5 + 1));

jac(3, 3) = (K5*g5*m4)/((km6 + m4)^2*(F11*m5 + 1)) - (K8*g8*m31*(F20*m4 + 1))/(km9 + m4*m31) - (K5*g5)/((km6 + m4)*(F11*m5 + 1)) + (K8*g8*m4*m31^2*(F20*m4 + 1))/(km9 + m4*m31)^2 - (K6*g6*m31*(F15*m6 + 1))/((km7 + m4*m31)*(F16*m31 + 1)) - (F23*K9*g9*m6)/((km10 + m6)*(F23*m4 + 1)^2) - (F20*K8*g8*m4*m31)/(km9 + m4*m31) + (K6*g6*m4*m31^2*(F15*m6 + 1))/((km7 + m4*m31)^2*(F16*m31 + 1));

jac(3, 4) = (K7*g7)/((km8 + m5)*(F19*m6 + 1)) - (K7*g7*m5)/((km8 + m5)^2*(F19*m6 + 1)) - (F8*K5*g5*m2)/((km5 + m2)*(F8*m5 + 1)^2) + (F11*K5*g5*m4)/((km6 + m4)*(F11*m5 + 1)^2);

jac(3, 5) = (K9*g9)/((km10 + m6)*(F23*m4 + 1)) - (K9*g9*m6)/((km10 + m6)^2*(F23*m4 + 1)) - (F19*K7*g7*m5)/((km8 + m5)*(F19*m6 + 1)^2) - (F15*K6*g6*m4*m31)/((km7 + m4*m31)*(F16*m31 + 1));

jac(3, 14) = (K8*g8*m4^2*m31*(F20*m4 + 1))/(km9 + m4*m31)^2 - (K8*g8*m4*(F20*m4 + 1))/(km9 + m4*m31) - (K6*g6*m4*(F15*m6 + 1))/((km7 + m4*m31)*(F16*m31 + 1)) + (K6*g6*m4^2*m31*(F15*m6 + 1))/((km7 + m4*m31)^2*(F16*m31 + 1)) + (F16*K6*g6*m4*m31*(F15*m6 + 1))/((km7 + m4*m31)*(F16*m31 + 1)^2);

jac(3, 20) = (K5*m2)/((km5 + m2)*(F8*m5 + 1)) - (K5*m4)/((km6 + m4)*(F11*m5 + 1));

jac(3, 21) = -(K6*m4*m31*(F15*m6 + 1))/((km7 + m4*m31)*(F16*m31 + 1));

jac(3, 22) = (K7*m5)/((km8 + m5)*(F19*m6 + 1));

jac(3, 23) = -(K8*m4*m31*(F20*m4 + 1))/(km9 + m4*m31);

jac(3, 24) = (K9*m6)/((km10 + m6)*(F23*m4 + 1));

jac(3, 35) = -(K5*g5*m2)/((km5 + m2)^2*(F8*m5 + 1));

jac(3, 36) = (K5*g5*m4)/((km6 + m4)^2*(F11*m5 + 1));

jac(3, 37) = (K6*g6*m4*m31*(F15*m6 + 1))/((km7 + m4*m31)^2*(F16*m31 + 1));

jac(3, 38) = -(K7*g7*m5)/((km8 + m5)^2*(F19*m6 + 1));

jac(3, 39) = (K8*g8*m4*m31*(F20*m4 + 1))/(km9 + m4*m31)^2;

jac(3, 40) = -(K9*g9*m6)/((km10 + m6)^2*(F23*m4 + 1));

jac(3, 57) = (g5*m2)/((km5 + m2)*(F8*m5 + 1)) - (g5*m4)/((km6 + m4)*(F11*m5 + 1));

jac(3, 58) = -(g6*m4*m31*(F15*m6 + 1))/((km7 + m4*m31)*(F16*m31 + 1));

jac(3, 59) = -(g8*m4*m31*(F20*m4 + 1))/(km9 + m4*m31);

jac(3, 71) = -(K5*g5*m2*m5)/((km5 + m2)*(F8*m5 + 1)^2);

jac(3, 72) = (K5*g5*m4*m5)/((km6 + m4)*(F11*m5 + 1)^2);

jac(3, 73) = -(K6*g6*m4*m6*m31)/((km7 + m4*m31)*(F16*m31 + 1));

jac(3, 74) = (K6*g6*m4*m31^2*(F15*m6 + 1))/((km7 + m4*m31)*(F16*m31 + 1)^2);

jac(3, 75) = -(K7*g7*m5*m6)/((km8 + m5)*(F19*m6 + 1)^2);

jac(3, 76) = -(K8*g8*m4^2*m31)/(km9 + m4*m31);

jac(3, 77) = -(K9*g9*m4*m6)/((km10 + m6)*(F23*m4 + 1)^2);

jac(4, 3) = (K6*g6*km7*m31*(F15*m6 + 1))/((km7 + m4*m31)^2*(F16*m31 + 1));

jac(4, 4) = (K10*g10*m5)/(km11 + m5)^2 - (K10*g10)/(km11 + m5) - (K7*g7)/((km8 + m5)*(F19*m6 + 1)) + (K7*g7*m5)/((km8 + m5)^2*(F19*m6 + 1));

jac(4, 5) = (F19*K7*g7*m5)/((km8 + m5)*(F19*m6 + 1)^2) + (F15*K6*g6*m4*m31)/((km7 + m4*m31)*(F16*m31 + 1));

jac(4, 6) = (K10*g10*km12*m8)/(km12 + m7*m8)^2;

jac(4, 7) = (K10*g10*km12*m7)/(km12 + m7*m8)^2;

jac(4, 14) = (K6*g6*m4*(F15*m6 + 1)*(km7 - F16*m4*m31^2))/((km7 + m4*m31)^2*(F16*m31 + 1)^2);

jac(4, 21) = (K6*m4*m31*(F15*m6 + 1))/((km7 + m4*m31)*(F16*m31 + 1));

jac(4, 22) = -(K7*m5)/((km8 + m5)*(F19*m6 + 1));

jac(4, 25) = (K10*m7*m8)/(km12 + m7*m8) - (K10*m5)/(km11 + m5);

jac(4, 37) = -(K6*g6*m4*m31*(F15*m6 + 1))/((km7 + m4*m31)^2*(F16*m31 + 1));

jac(4, 38) = (K7*g7*m5)/((km8 + m5)^2*(F19*m6 + 1));

jac(4, 41) = (K10*g10*m5)/(km11 + m5)^2;

jac(4, 42) = -(K10*g10*m7*m8)/(km12 + m7*m8)^2;

jac(4, 58) = (g6*m4*m31*(F15*m6 + 1))/((km7 + m4*m31)*(F16*m31 + 1));

jac(4, 60) = (g10*m7*m8)/(km12 + m7*m8) - (g10*m5)/(km11 + m5);

jac(4, 73) = (K6*g6*m4*m6*m31)/((km7 + m4*m31)*(F16*m31 + 1));

jac(4, 74) = -(K6*g6*m4*m31^2*(F15*m6 + 1))/((km7 + m4*m31)*(F16*m31 + 1)^2);

jac(4, 75) = (K7*g7*m5*m6)/((km8 + m5)*(F19*m6 + 1)^2);

jac(5, 3) = (K8*g8*m31*(F20*m4 + 1))/(km9 + m4*m31) - (K8*g8*m4*m31^2*(F20*m4 + 1))/(km9 + m4*m31)^2 + (F23*K9*g9*m6)/((km10 + m6)*(F23*m4 + 1)^2) + (F20*K8*g8*m4*m31)/(km9 + m4*m31);

jac(5, 5) = -(K9*g9*km10)/((km10 + m6)^2*(F23*m4 + 1));

jac(5, 14) = (K8*g8*km9*m4*(F20*m4 + 1))/(km9 + m4*m31)^2;

jac(5, 23) = (K8*m4*m31*(F20*m4 + 1))/(km9 + m4*m31);

jac(5, 24) = -(K9*m6)/((km10 + m6)*(F23*m4 + 1));

jac(5, 39) = -(K8*g8*m4*m31*(F20*m4 + 1))/(km9 + m4*m31)^2;

jac(5, 40) = (K9*g9*m6)/((km10 + m6)^2*(F23*m4 + 1));

jac(5, 59) = (g8*m4*m31*(F20*m4 + 1))/(km9 + m4*m31);

jac(5, 76) = (K8*g8*m4^2*m31)/(km9 + m4*m31);

jac(5, 77) = (K9*g9*m4*m6)/((km10 + m6)*(F23*m4 + 1)^2);

jac(6, 4) = (K10*g10*km11)/(km11 + m5)^2;

jac(6, 6) = (K11*g11*m7)/(km13 + m7)^2 - (K11*g11)/(km13 + m7) - (K10*g10*m8)/(km12 + m7*m8) + (K10*g10*m7*m8^2)/(km12 + m7*m8)^2;

jac(6, 7) = (K11*g11)/(km14 + m8) - (K11*g11*m8)/(km14 + m8)^2 - (K10*g10*m7)/(km12 + m7*m8) + (K10*g10*m7^2*m8)/(km12 + m7*m8)^2;

jac(6, 25) = (K10*m5)/(km11 + m5) - (K10*m7*m8)/(km12 + m7*m8);

jac(6, 26) = (K11*m8)/(km14 + m8) - (K11*m7)/(km13 + m7);

jac(6, 41) = -(K10*g10*m5)/(km11 + m5)^2;

jac(6, 42) = (K10*g10*m7*m8)/(km12 + m7*m8)^2;

jac(6, 43) = (K11*g11*m7)/(km13 + m7)^2;

jac(6, 44) = -(K11*g11*m8)/(km14 + m8)^2;

jac(6, 60) = (g10*m5)/(km11 + m5) - (g10*m7*m8)/(km12 + m7*m8);

jac(6, 61) = (g11*m8)/(km14 + m8) - (g11*m7)/(km13 + m7);

jac(7, 4) = (K10*g10*km11)/(km11 + m5)^2;

jac(7, 6) = (K11*g11)/(km13 + m7) - (K11*g11*m7)/(km13 + m7)^2 - (K10*g10*m8)/(km12 + m7*m8) + (K10*g10*m7*m8^2)/(km12 + m7*m8)^2;

jac(7, 7) = (K11*g11*m8)/(km14 + m8)^2 - (K11*g11)/(km14 + m8) - (K10*g10*m7)/(km12 + m7*m8) - (K12*g12*m34)/(km15 + m8*m34) + (K10*g10*m7^2*m8)/(km12 + m7*m8)^2 + (K12*g12*m8*m34^2)/(km15 + m8*m34)^2;

jac(7, 8) = (K12*g12*km16*m33)/(km16 + m9*m33)^2;

jac(7, 16) = (K12*g12*km16*m9)/(km16 + m9*m33)^2;

jac(7, 17) = -(K12*g12*km15*m8)/(km15 + m8*m34)^2;

jac(7, 25) = (K10*m5)/(km11 + m5) - (K10*m7*m8)/(km12 + m7*m8);

jac(7, 26) = (K11*m7)/(km13 + m7) - (K11*m8)/(km14 + m8);

jac(7, 27) = (K12*m9*m33)/(km16 + m9*m33) - (K12*m8*m34)/(km15 + m8*m34);

jac(7, 41) = -(K10*g10*m5)/(km11 + m5)^2;

jac(7, 42) = (K10*g10*m7*m8)/(km12 + m7*m8)^2;

jac(7, 43) = -(K11*g11*m7)/(km13 + m7)^2;

jac(7, 44) = (K11*g11*m8)/(km14 + m8)^2;

jac(7, 45) = (K12*g12*m8*m34)/(km15 + m8*m34)^2;

jac(7, 46) = -(K12*g12*m9*m33)/(km16 + m9*m33)^2;

jac(7, 60) = (g10*m5)/(km11 + m5) - (g10*m7*m8)/(km12 + m7*m8);

jac(7, 61) = (g11*m7)/(km13 + m7) - (g11*m8)/(km14 + m8);

jac(7, 62) = (g12*m9*m33)/(km16 + m9*m33) - (g12*m8*m34)/(km15 + m8*m34);

jac(8, 7) = (K12*g12*km15*m34)/(km15 + m8*m34)^2;

jac(8, 8) = (K12*g12*m9*m33^2)/(km16 + m9*m33)^2 - (K13*g13*m32)/(km17 + m9*m32) - (K12*g12*m33)/(km16 + m9*m33) + (K13*g13*m9*m32^2)/(km17 + m9*m32)^2;

jac(8, 9) = (K13*g13*km18)/(km18 + m10)^2;

jac(8, 15) = -(K13*g13*km17*m9)/(km17 + m9*m32)^2;

jac(8, 16) = -(K12*g12*km16*m9)/(km16 + m9*m33)^2;

jac(8, 17) = (K12*g12*km15*m8)/(km15 + m8*m34)^2;

jac(8, 27) = (K12*m8*m34)/(km15 + m8*m34) - (K12*m9*m33)/(km16 + m9*m33);

jac(8, 28) = (K13*m10)/(km18 + m10) - (K13*m9*m32)/(km17 + m9*m32);

jac(8, 45) = -(K12*g12*m8*m34)/(km15 + m8*m34)^2;

jac(8, 46) = (K12*g12*m9*m33)/(km16 + m9*m33)^2;

jac(8, 47) = (K13*g13*m9*m32)/(km17 + m9*m32)^2;

jac(8, 48) = -(K13*g13*m10)/(km18 + m10)^2;

jac(8, 62) = (g12*m8*m34)/(km15 + m8*m34) - (g12*m9*m33)/(km16 + m9*m33);

jac(8, 63) = (g13*m10)/(km18 + m10) - (g13*m9*m32)/(km17 + m9*m32);

jac(9, 8) = (K13*g13*km17*m32)/(km17 + m9*m32)^2;

jac(9, 9) = (K13*g13*m10)/(km18 + m10)^2 - (K14*g14)/(km19 + m10) - (K13*g13)/(km18 + m10) + (K14*g14*m10)/(km19 + m10)^2;

jac(9, 10) = (K14*g14*km20)/(km20 + m11)^2;

jac(9, 15) = (K13*g13*km17*m9)/(km17 + m9*m32)^2;

jac(9, 28) = (K13*m9*m32)/(km17 + m9*m32) - (K13*m10)/(km18 + m10);

jac(9, 29) = (K14*m11)/(km20 + m11) - (K14*m10)/(km19 + m10);

jac(9, 47) = -(K13*g13*m9*m32)/(km17 + m9*m32)^2;

jac(9, 48) = (K13*g13*m10)/(km18 + m10)^2;

jac(9, 49) = (K14*g14*m10)/(km19 + m10)^2;

jac(9, 50) = -(K14*g14*m11)/(km20 + m11)^2;

jac(9, 63) = (g13*m9*m32)/(km17 + m9*m32) - (g13*m10)/(km18 + m10);

jac(9, 64) = (g14*m11)/(km20 + m11) - (g14*m10)/(km19 + m10);

jac(10, 9) = (K14*g14*km19)/(km19 + m10)^2;

jac(10, 10) = (K14*g14*m11)/(km20 + m11)^2 - (K15*g15)/(km21 + m11) - (K14*g14)/(km20 + m11) + (K15*g15*m11)/(km21 + m11)^2;

jac(10, 11) = (K15*g15*km22)/(km22 + m12)^2;

jac(10, 29) = (K14*m10)/(km19 + m10) - (K14*m11)/(km20 + m11);

jac(10, 30) = (K15*m12)/(km22 + m12) - (K15*m11)/(km21 + m11);

jac(10, 49) = -(K14*g14*m10)/(km19 + m10)^2;

jac(10, 50) = (K14*g14*m11)/(km20 + m11)^2;

jac(10, 51) = (K15*g15*m11)/(km21 + m11)^2;

jac(10, 52) = -(K15*g15*m12)/(km22 + m12)^2;

jac(10, 64) = (g14*m10)/(km19 + m10) - (g14*m11)/(km20 + m11);

jac(10, 65) = (g15*m12)/(km22 + m12) - (g15*m11)/(km21 + m11);

jac(11, 4) = -(F27*K16*g16*m12*m32)/((km23 + m12*m32)*(F28*m31 + 1));

jac(11, 10) = (K15*g15*km21)/(km21 + m11)^2;

jac(11, 11) = (K15*g15*m12)/(km22 + m12)^2 - (K15*g15)/(km22 + m12) - (K16*g16*m32*(F27*m5 + 1))/((km23 + m12*m32)*(F28*m31 + 1)) + (K16*g16*m12*m32^2*(F27*m5 + 1))/((km23 + m12*m32)^2*(F28*m31 + 1));

jac(11, 14) = (F28*K16*g16*m12*m32*(F27*m5 + 1))/((km23 + m12*m32)*(F28*m31 + 1)^2);

jac(11, 15) = -(K16*g16*km23*m12*(F27*m5 + 1))/((km23 + m12*m32)^2*(F28*m31 + 1));

jac(11, 30) = (K15*m11)/(km21 + m11) - (K15*m12)/(km22 + m12);

jac(11, 31) = -(K16*m12*m32*(F27*m5 + 1))/((km23 + m12*m32)*(F28*m31 + 1));

jac(11, 51) = -(K15*g15*m11)/(km21 + m11)^2;

jac(11, 52) = (K15*g15*m12)/(km22 + m12)^2;

jac(11, 53) = (K16*g16*m12*m32*(F27*m5 + 1))/((km23 + m12*m32)^2*(F28*m31 + 1));

jac(11, 65) = (g15*m11)/(km21 + m11) - (g15*m12)/(km22 + m12);

jac(11, 66) = -(g16*m12*m32*(F27*m5 + 1))/((km23 + m12*m32)*(F28*m31 + 1));

jac(11, 78) = -(K16*g16*m5*m12*m32)/((km23 + m12*m32)*(F28*m31 + 1));

jac(11, 79) = (K16*g16*m12*m31*m32*(F27*m5 + 1))/((km23 + m12*m32)*(F28*m31 + 1)^2);

jac(12, 4) = (F27*K16*g16*m12*m32)/((km23 + m12*m32)*(F28*m31 + 1));

jac(12, 11) = (K16*g16*km23*m32*(F27*m5 + 1))/((km23 + m12*m32)^2*(F28*m31 + 1));

jac(12, 12) = -(K17*g17*km24*m33)/(km24 + m13*m33)^2;

jac(12, 14) = -(F28*K16*g16*m12*m32*(F27*m5 + 1))/((km23 + m12*m32)*(F28*m31 + 1)^2);

jac(12, 15) = (K16*g16*km23*m12*(F27*m5 + 1))/((km23 + m12*m32)^2*(F28*m31 + 1));

jac(12, 16) = -(K17*g17*km24*m13)/(km24 + m13*m33)^2;

jac(12, 31) = (K16*m12*m32*(F27*m5 + 1))/((km23 + m12*m32)*(F28*m31 + 1));

jac(12, 32) = -(K17*m13*m33)/(km24 + m13*m33);

jac(12, 53) = -(K16*g16*m12*m32*(F27*m5 + 1))/((km23 + m12*m32)^2*(F28*m31 + 1));

jac(12, 54) = (K17*g17*m13*m33)/(km24 + m13*m33)^2;

jac(12, 66) = (g16*m12*m32*(F27*m5 + 1))/((km23 + m12*m32)*(F28*m31 + 1));

jac(12, 67) = -(g17*m13*m33)/(km24 + m13*m33);

jac(12, 78) = (K16*g16*m5*m12*m32)/((km23 + m12*m32)*(F28*m31 + 1));

jac(12, 79) = -(K16*g16*m12*m31*m32*(F27*m5 + 1))/((km23 + m12*m32)*(F28*m31 + 1)^2);

jac(13, 12) = (K17*g17*km24*m33)/(km24 + m13*m33)^2;

jac(13, 16) = (K17*g17*km24*m13)/(km24 + m13*m33)^2;

jac(13, 32) = (K17*m13*m33)/(km24 + m13*m33);

jac(13, 54) = -(K17*g17*m13*m33)/(km24 + m13*m33)^2;

jac(13, 67) = (g17*m13*m33)/(km24 + m13*m33);

jac(14, 1) = (F5*K3*g3*m3*m31*(F4*m32 + 1))/((km3 + m3*m31)*(F5*m2 + 1)^2);

jac(14, 2) = -(K3*g3*km3*m31*(F4*m32 + 1))/((km3 + m3*m31)^2*(F5*m2 + 1));

jac(14, 3) = (K8*g8*m4*m31^2*(F20*m4 + 1))/(km9 + m4*m31)^2 - (K8*g8*m31*(F20*m4 + 1))/(km9 + m4*m31) - (K6*g6*m31*(F15*m6 + 1))/((km7 + m4*m31)*(F16*m31 + 1)) - (F20*K8*g8*m4*m31)/(km9 + m4*m31) + (K6*g6*m4*m31^2*(F15*m6 + 1))/((km7 + m4*m31)^2*(F16*m31 + 1));

jac(14, 4) = (F27*K16*g16*m12*m32)/((km23 + m12*m32)*(F28*m31 + 1));

jac(14, 5) = -(F15*K6*g6*m4*m31)/((km7 + m4*m31)*(F16*m31 + 1));

jac(14, 8) = (K13*g13*km17*m32)/(km17 + m9*m32)^2;

jac(14, 11) = (K16*g16*km23*m32*(F27*m5 + 1))/((km23 + m12*m32)^2*(F28*m31 + 1));

jac(14, 14) = (K8*g8*m4^2*m31*(F20*m4 + 1))/(km9 + m4*m31)^2 - (K8*g8*m4*(F20*m4 + 1))/(km9 + m4*m31) - (K3*g3*m3*(F4*m32 + 1))/((km3 + m3*m31)*(F5*m2 + 1)) - (K6*g6*m4*(F15*m6 + 1))/((km7 + m4*m31)*(F16*m31 + 1)) + (K3*g3*m3^2*m31*(F4*m32 + 1))/((km3 + m3*m31)^2*(F5*m2 + 1)) + (K6*g6*m4^2*m31*(F15*m6 + 1))/((km7 + m4*m31)^2*(F16*m31 + 1)) + (F16*K6*g6*m4*m31*(F15*m6 + 1))/((km7 + m4*m31)*(F16*m31 + 1)^2) - (F28*K16*g16*m12*m32*(F27*m5 + 1))/((km23 + m12*m32)*(F28*m31 + 1)^2);

jac(14, 15) = (K13*g13*m9)/(km17 + m9*m32) - (K13*g13*m9^2*m32)/(km17 + m9*m32)^2 + (K16*g16*m12*(F27*m5 + 1))/((km23 + m12*m32)*(F28*m31 + 1)) - (F4*K3*g3*m3*m31)/((km3 + m3*m31)*(F5*m2 + 1)) - (K16*g16*m12^2*m32*(F27*m5 + 1))/((km23 + m12*m32)^2*(F28*m31 + 1));

jac(14, 18) = -(K3*m3*m31*(F4*m32 + 1))/((km3 + m3*m31)*(F5*m2 + 1));

jac(14, 21) = -(K6*m4*m31*(F15*m6 + 1))/((km7 + m4*m31)*(F16*m31 + 1));

jac(14, 23) = -(K8*m4*m31*(F20*m4 + 1))/(km9 + m4*m31);

jac(14, 28) = (K13*m9*m32)/(km17 + m9*m32);

jac(14, 31) = (K16*m12*m32*(F27*m5 + 1))/((km23 + m12*m32)*(F28*m31 + 1));

jac(14, 33) = (K3*g3*m3*m31*(F4*m32 + 1))/((km3 + m3*m31)^2*(F5*m2 + 1));

jac(14, 37) = (K6*g6*m4*m31*(F15*m6 + 1))/((km7 + m4*m31)^2*(F16*m31 + 1));

jac(14, 39) = (K8*g8*m4*m31*(F20*m4 + 1))/(km9 + m4*m31)^2;

jac(14, 47) = -(K13*g13*m9*m32)/(km17 + m9*m32)^2;

jac(14, 53) = -(K16*g16*m12*m32*(F27*m5 + 1))/((km23 + m12*m32)^2*(F28*m31 + 1));

jac(14, 55) = -(g3*m3*m31*(F4*m32 + 1))/((km3 + m3*m31)*(F5*m2 + 1));

jac(14, 58) = -(g6*m4*m31*(F15*m6 + 1))/((km7 + m4*m31)*(F16*m31 + 1));

jac(14, 59) = -(g8*m4*m31*(F20*m4 + 1))/(km9 + m4*m31);

jac(14, 63) = (g13*m9*m32)/(km17 + m9*m32);

jac(14, 66) = (g16*m12*m32*(F27*m5 + 1))/((km23 + m12*m32)*(F28*m31 + 1));

jac(14, 68) = -(K3*g3*m3*m31*m32)/((km3 + m3*m31)*(F5*m2 + 1));

jac(14, 69) = (K3*g3*m2*m3*m31*(F4*m32 + 1))/((km3 + m3*m31)*(F5*m2 + 1)^2);

jac(14, 73) = -(K6*g6*m4*m6*m31)/((km7 + m4*m31)*(F16*m31 + 1));

jac(14, 74) = (K6*g6*m4*m31^2*(F15*m6 + 1))/((km7 + m4*m31)*(F16*m31 + 1)^2);

jac(14, 76) = -(K8*g8*m4^2*m31)/(km9 + m4*m31);

jac(14, 78) = (K16*g16*m5*m12*m32)/((km23 + m12*m32)*(F28*m31 + 1));

jac(14, 79) = -(K16*g16*m12*m31*m32*(F27*m5 + 1))/((km23 + m12*m32)*(F28*m31 + 1)^2);

jac(15, 1) = -(F5*K3*g3*m3*m31*(F4*m32 + 1))/((km3 + m3*m31)*(F5*m2 + 1)^2);

jac(15, 2) = (K3*g3*km3*m31*(F4*m32 + 1))/((km3 + m3*m31)^2*(F5*m2 + 1));

jac(15, 3) = (K8*g8*m31*(F20*m4 + 1))/(km9 + m4*m31) - (K8*g8*m4*m31^2*(F20*m4 + 1))/(km9 + m4*m31)^2 + (K6*g6*m31*(F15*m6 + 1))/((km7 + m4*m31)*(F16*m31 + 1)) + (F20*K8*g8*m4*m31)/(km9 + m4*m31) - (K6*g6*m4*m31^2*(F15*m6 + 1))/((km7 + m4*m31)^2*(F16*m31 + 1));

jac(15, 4) = -(F27*K16*g16*m12*m32)/((km23 + m12*m32)*(F28*m31 + 1));

jac(15, 5) = (F15*K6*g6*m4*m31)/((km7 + m4*m31)*(F16*m31 + 1));

jac(15, 8) = -(K13*g13*km17*m32)/(km17 + m9*m32)^2;

jac(15, 11) = -(K16*g16*km23*m32*(F27*m5 + 1))/((km23 + m12*m32)^2*(F28*m31 + 1));

jac(15, 14) = (K8*g8*m4*(F20*m4 + 1))/(km9 + m4*m31) - (K8*g8*m4^2*m31*(F20*m4 + 1))/(km9 + m4*m31)^2 + (K3*g3*m3*(F4*m32 + 1))/((km3 + m3*m31)*(F5*m2 + 1)) + (K6*g6*m4*(F15*m6 + 1))/((km7 + m4*m31)*(F16*m31 + 1)) - (K3*g3*m3^2*m31*(F4*m32 + 1))/((km3 + m3*m31)^2*(F5*m2 + 1)) - (K6*g6*m4^2*m31*(F15*m6 + 1))/((km7 + m4*m31)^2*(F16*m31 + 1)) - (F16*K6*g6*m4*m31*(F15*m6 + 1))/((km7 + m4*m31)*(F16*m31 + 1)^2) + (F28*K16*g16*m12*m32*(F27*m5 + 1))/((km23 + m12*m32)*(F28*m31 + 1)^2);

jac(15, 15) = (K13*g13*m9^2*m32)/(km17 + m9*m32)^2 - (K13*g13*m9)/(km17 + m9*m32) - (K16*g16*m12*(F27*m5 + 1))/((km23 + m12*m32)*(F28*m31 + 1)) + (F4*K3*g3*m3*m31)/((km3 + m3*m31)*(F5*m2 + 1)) + (K16*g16*m12^2*m32*(F27*m5 + 1))/((km23 + m12*m32)^2*(F28*m31 + 1));

jac(15, 18) = (K3*m3*m31*(F4*m32 + 1))/((km3 + m3*m31)*(F5*m2 + 1));

jac(15, 21) = (K6*m4*m31*(F15*m6 + 1))/((km7 + m4*m31)*(F16*m31 + 1));

jac(15, 23) = (K8*m4*m31*(F20*m4 + 1))/(km9 + m4*m31);

jac(15, 28) = -(K13*m9*m32)/(km17 + m9*m32);

jac(15, 31) = -(K16*m12*m32*(F27*m5 + 1))/((km23 + m12*m32)*(F28*m31 + 1));

jac(15, 33) = -(K3*g3*m3*m31*(F4*m32 + 1))/((km3 + m3*m31)^2*(F5*m2 + 1));

jac(15, 37) = -(K6*g6*m4*m31*(F15*m6 + 1))/((km7 + m4*m31)^2*(F16*m31 + 1));

jac(15, 39) = -(K8*g8*m4*m31*(F20*m4 + 1))/(km9 + m4*m31)^2;

jac(15, 47) = (K13*g13*m9*m32)/(km17 + m9*m32)^2;

jac(15, 53) = (K16*g16*m12*m32*(F27*m5 + 1))/((km23 + m12*m32)^2*(F28*m31 + 1));

jac(15, 55) = (g3*m3*m31*(F4*m32 + 1))/((km3 + m3*m31)*(F5*m2 + 1));

jac(15, 58) = (g6*m4*m31*(F15*m6 + 1))/((km7 + m4*m31)*(F16*m31 + 1));

jac(15, 59) = (g8*m4*m31*(F20*m4 + 1))/(km9 + m4*m31);

jac(15, 63) = -(g13*m9*m32)/(km17 + m9*m32);

jac(15, 66) = -(g16*m12*m32*(F27*m5 + 1))/((km23 + m12*m32)*(F28*m31 + 1));

jac(15, 68) = (K3*g3*m3*m31*m32)/((km3 + m3*m31)*(F5*m2 + 1));

jac(15, 69) = -(K3*g3*m2*m3*m31*(F4*m32 + 1))/((km3 + m3*m31)*(F5*m2 + 1)^2);

jac(15, 73) = (K6*g6*m4*m6*m31)/((km7 + m4*m31)*(F16*m31 + 1));

jac(15, 74) = -(K6*g6*m4*m31^2*(F15*m6 + 1))/((km7 + m4*m31)*(F16*m31 + 1)^2);

jac(15, 76) = (K8*g8*m4^2*m31)/(km9 + m4*m31);

jac(15, 78) = -(K16*g16*m5*m12*m32)/((km23 + m12*m32)*(F28*m31 + 1));

jac(15, 79) = (K16*g16*m12*m31*m32*(F27*m5 + 1))/((km23 + m12*m32)*(F28*m31 + 1)^2);

jac(16, 7) = (K12*g12*km15*m34)/(km15 + m8*m34)^2;

jac(16, 8) = -(K12*g12*km16*m33)/(km16 + m9*m33)^2;

jac(16, 12) = -(K17*g17*km24*m33)/(km24 + m13*m33)^2;

jac(16, 16) = (K12*g12*m9^2*m33)/(km16 + m9*m33)^2 - (K17*g17*m13)/(km24 + m13*m33) - (K12*g12*m9)/(km16 + m9*m33) + (K17*g17*m13^2*m33)/(km24 + m13*m33)^2;

jac(16, 17) = (K12*g12*km15*m8)/(km15 + m8*m34)^2;

jac(16, 27) = (K12*m8*m34)/(km15 + m8*m34) - (K12*m9*m33)/(km16 + m9*m33);

jac(16, 32) = -(K17*m13*m33)/(km24 + m13*m33);

jac(16, 45) = -(K12*g12*m8*m34)/(km15 + m8*m34)^2;

jac(16, 46) = (K12*g12*m9*m33)/(km16 + m9*m33)^2;

jac(16, 54) = (K17*g17*m13*m33)/(km24 + m13*m33)^2;

jac(16, 62) = (g12*m8*m34)/(km15 + m8*m34) - (g12*m9*m33)/(km16 + m9*m33);

jac(16, 67) = -(g17*m13*m33)/(km24 + m13*m33);

jac(17, 7) = -(K12*g12*km15*m34)/(km15 + m8*m34)^2;

jac(17, 8) = (K12*g12*km16*m33)/(km16 + m9*m33)^2;

jac(17, 12) = (K17*g17*km24*m33)/(km24 + m13*m33)^2;

jac(17, 16) = (K12*g12*m9)/(km16 + m9*m33) + (K17*g17*m13)/(km24 + m13*m33) - (K12*g12*m9^2*m33)/(km16 + m9*m33)^2 - (K17*g17*m13^2*m33)/(km24 + m13*m33)^2;

jac(17, 17) = -(K12*g12*km15*m8)/(km15 + m8*m34)^2;

jac(17, 27) = (K12*m9*m33)/(km16 + m9*m33) - (K12*m8*m34)/(km15 + m8*m34);

jac(17, 32) = (K17*m13*m33)/(km24 + m13*m33);

jac(17, 45) = (K12*g12*m8*m34)/(km15 + m8*m34)^2;

jac(17, 46) = -(K12*g12*m9*m33)/(km16 + m9*m33)^2;

jac(17, 54) = -(K17*g17*m13*m33)/(km24 + m13*m33)^2;

jac(17, 62) = (g12*m9*m33)/(km16 + m9*m33) - (g12*m8*m34)/(km15 + m8*m34);

jac(17, 67) = (g17*m13*m33)/(km24 + m13*m33);

jac(18, 18) = 1;

jac(19, 19) = 1;

jac(20, 20) = 1;

jac(21, 21) = 1;

jac(22, 22) = 1;

jac(23, 23) = 1;

jac(24, 24) = 1;

jac(25, 25) = 1;

jac(26, 26) = 1;

jac(27, 27) = 1;

jac(28, 28) = 1;

jac(29, 29) = 1;

jac(30, 30) = 1;

jac(31, 31) = 1;

jac(32, 32) = 1;



