function [xdot] = Glycolysis_extended_model( t,x_v)

NOM = 79;              % the number of molecules
xdot = zeros(NOM, 1);


%%%%% METABOLITES INITIAL CONCENTRATION in nM %%%%%
%%% m39 is not included intentionally


m2	=	x_v(1);     %	glc_6P					
m3	=	x_v(2);     %	glc				
m4	=	x_v(3);     %	f_6P					
m5	=	x_v(4);     %	f_16_BP					
m6	=	x_v(5);     %	f_26_BP				
m7	=	x_v(6);     %	dhap				
m8	=	x_v(7);     %	ga_3P					
m9	=	x_v(8);     %	BPG_13					
m10	=	x_v(9);	%	PG_3					
m11	=	x_v(10);	%	PG_2					
m12	=	x_v(11);	%	PEP					
m13	=	x_v(12);	%	PYR				
m14	=	x_v(13);	%	LACTATE		
				
m31	=	x_v(14);	%	ATP					
m32	=	x_v(15);	%	ADP			
m33	=	x_v(16);	%	NADH					
m34	=	x_v(17);	%	NAD	



%%%%% ENZYME/GENE INITIAL CONCENTRATION g in nM  %%%%%

g3	=	x_v(18);	%Hexokinase				
g4	=	x_v(19);	%Glucose-6-phosphatase				
g5	=	x_v(20);	%Phosphoglucoisomerase				
g6	=	x_v(21);	%Phosphofructokinase_1				
g7	=	x_v(22);	%Fructose-1,6-bisphosphatase				
g8	=	x_v(23);	%Phosphofructokinase_2				
g9	=	x_v(24);	%Fructose-2,6-bisphosphatase				
g10	=	x_v(25);	%Aldolase
g11	=	x_v(26);	%Triose_phosphate_isomerase				
g12	=	x_v(27);	%Glyceraldehyde-3-phosphate_dehydrogenase				
g13	=	x_v(28);	%Phosphoglycerate_kinase				
g14	=	x_v(29);	%Phosphoglycerate_mutase				
g15	=	x_v(30);	%Enolase				
g16	=	x_v(31);	%Pyruvate_kinase				
g17	=	x_v(32);	%Lactate_dehydrogenase	


% g3 2,3BPG	


km3	 =	x_v(33);
km4	 =	x_v(34);
km5	 =	x_v(35);
km6	 =	x_v(36);
km7	 =	x_v(37);
km8	 =	x_v(38);
km9	 =	x_v(39);
km10	 =	x_v(40);
km11	 =	x_v(41);
km12	 =	x_v(42);
km13	 =	x_v(43);
km14	 =	x_v(44);
km15	 =	x_v(45);
km16	 =	x_v(46);
km17	 =	x_v(47);
km18	 =	x_v(48);
km19	 =	x_v(49);
km20	 =	x_v(50);
km21	 =	x_v(51);
km22	 =	x_v(52);
km23	 =	x_v(53);
km24	 =	x_v(54);


%%%%% METABOLIC REACTION ENZYME RATE CONATANT K in S^-1  %%%%%


K3	 =	x_v(55);	% hk
K4	 =	x_v(56);	% g6Pase
K5	 =	x_v(57);	% pgi
K6	 =	x_v(58);	% pfk1
K7	 =	0;%0.017	;	% f16Bpase
K8	 =	x_v(59);	% pfk2
K9	 =	0;%0.0766	;	% f26Bpase
K10	 =	x_v(60);	% ald
K11	 =	x_v(61);	% tpi
K12	 =	x_v(62);	% gcld3PD
K13	 =	x_v(63);	% pglc_kn
K14	 =	x_v(64);	% pglc_m
K15	 =	x_v(65);	% enl
K16	 =	x_v(66);	% pyrk
K17	 =	x_v(67);	% lacd


F4	=	x_v(68);
F5	=	x_v(69);
F7	=	x_v(70);
F8	=	x_v(71);
F11	=	x_v(72);
F15	=	x_v(73);
F16	=	x_v(74);
F19	=	x_v(75);
F20	=	x_v(76);
F23	=	x_v(77);
F27	=	x_v(78);
F28	=	x_v(79);



%%%%% METABOLIC ODE PART


xdot(1) = - (K4*g4*m2*(F7*m2 + 1))/(km4 + m2)  - (K5*g5*m2)/((km5 + m2)*(F8*m5 + 1)) + (K5*g5*m4)/((km6 + m4)*(F11*m5 + 1)) + (K3*g3*m3*m31*(F4*m32 + 1))/((km3 + m3*m31)*(F5*m2 + 1));
xdot(2) =  (K4*g4*m2*(F7*m2 + 1))/(km4 + m2) - (K3*g3*m3*m31*(F4*m32 + 1))/((km3 + m3*m31)*(F5*m2 + 1));
xdot(3) =  (K9*g9*m6)/((km10 + m6)*(F23*m4 + 1)) + (K7*g7*m5)/((km8 + m5)*(F19*m6 + 1)) + (K5*g5*m2)/((km5 + m2)*(F8*m5 + 1)) - (K5*g5*m4)/((km6 + m4)*(F11*m5 + 1)) - (K8*g8*m4*m31*(F20*m4 + 1))/((km9 + m4*m31)) - (K6*g6*m4*m31*(F15*m6 + 1))/((km7 + m4*m31)*(F16*m31 + 1));
xdot(4) = (K10*g10*m7*m8)/(km12 + m7*m8) - (K10*g10*m5)/(km11 + m5) - (K7*g7*m5)/((km8 + m5)*(F19*m6 + 1)) + (K6*g6*m4*m31*(F15*m6 + 1))/((km7 + m4*m31)*(F16*m31 + 1));
xdot(5) = (K8*g8*m4*m31*(F20*m4 + 1))/((km9 + m4*m31)) - (K9*g9*m6)/((km10 + m6)*(F23*m4 + 1));
xdot(6) = (K10*g10*m5)/(km11 + m5) - (K11*g11*m7)/(km13 + m7) + (K11*g11*m8)/(km14 + m8) - (K10*g10*m7*m8)/(km12 + m7*m8);
xdot(7) = (K10*g10*m5)/(km11 + m5) + (K11*g11*m7)/(km13 + m7) - (K11*g11*m8)/(km14 + m8) - (K10*g10*m7*m8)/(km12 + m7*m8) - (K12*g12*m8*m34)/(km15 + m8*m34) + (K12*g12*m9*m33)/(km16 + m9*m33);
xdot(8) = (K13*g13*m10)/(km18 + m10) + (K12*g12*m8*m34)/(km15 + m8*m34) - (K12*g12*m9*m33)/(km16 + m9*m33) - (K13*g13*m9*m32)/(km17 + m9*m32);
xdot(9) = (K13*g13*m9*m32)/(km17 + m9*m32) - (K13*g13*m10)/(km18 + m10) - (K14*g14*m10)/(km19 + m10) + (K14*g14*m11)/(km20 + m11);
xdot(10) = (K15*g15*m12)/(km22 + m12) - (K15*g15*m11)/(km21 + m11) + (K14*g14*m10)/(km19 + m10) - (K14*g14*m11)/(km20 + m11);
xdot(11) = (K15*g15*m11)/(km21 + m11) - (K15*g15*m12)/(km22 + m12) - (K16*g16*m12*m32*(F27*m5 + 1))/((km23 + m12*m32)*(F28*m31 + 1));
xdot(12) = (K16*g16*m12*m32*(F27*m5 + 1))/((km23 + m12*m32)*(F28*m31 + 1)) - (K17*g17*m13*m33)/(km24 + m13*m33);
xdot(13) = (K17*g17*m13*m33)/(km24 + m13*m33);


xdot(14) = (K13*g13*m9*m32)/(km17 + m9*m32)  - (K8*g8*m4*m31*(F20*m4 + 1))/((km9 + m4*m31))  - (K3*g3*m3*m31*(F4*m32 + 1))/((km3 + m3*m31)*(F5*m2 + 1)) - (K6*g6*m4*m31*(F15*m6 + 1))/((km7 + m4*m31)*(F16*m31 + 1)) + (K16*g16*m12*m32*(F27*m5 + 1))/((km23 + m12*m32)*(F28*m31 + 1));
xdot(15) =  - (K13*g13*m9*m32)/(km17 + m9*m32) + (K8*g8*m4*m31*(F20*m4 + 1))/((km9 + m4*m31)) + (K3*g3*m3*m31*(F4*m32 + 1))/((km3 + m3*m31)*(F5*m2 + 1)) + (K6*g6*m4*m31*(F15*m6 + 1))/((km7 + m4*m31)*(F16*m31 + 1)) - (K16*g16*m12*m32*(F27*m5 + 1))/((km23 + m12*m32)*(F28*m31 + 1));
xdot(16) = (K12*g12*m8*m34)/(km15 + m8*m34) - (K12*g12*m9*m33)/(km16 + m9*m33) - (K17*g17*m13*m33)/(km24 + m13*m33);
xdot(17) = (K12*g12*m9*m33)/(km16 + m9*m33) - (K12*g12*m8*m34)/(km15 + m8*m34) + (K17*g17*m13*m33)/(km24 + m13*m33);
xdot(18) = g3;
xdot(19) = g4;
xdot(20) = g5;
xdot(21) = g6;
xdot(22) = g7;
xdot(23) = g8;
xdot(24) = g9;
xdot(25) = g10;
xdot(26) = g11;
xdot(27) = g12;
xdot(28) = g13;
xdot(29) = g14;
xdot(30) = g15;
xdot(31) = g16;
xdot(32) = g17;
xdot(33) = 0;
xdot(34) = 0;
xdot(35) = 0;
xdot(36) = 0;
xdot(37) = 0;
xdot(38) = 0;
xdot(39) = 0;
xdot(40) = 0;
xdot(41) = 0;
xdot(42) = 0;
xdot(43) = 0;
xdot(44) = 0;
xdot(45) = 0;
xdot(46) = 0;
xdot(47) = 0;
xdot(48) = 0;
xdot(49) = 0;
xdot(50) = 0;
xdot(51) = 0;
xdot(52) = 0;
xdot(53) = 0;
xdot(54) = 0;
xdot(55) = 0;
xdot(56) = 0;
xdot(57) = 0;
xdot(58) = 0;
xdot(59) = 0;
xdot(60) = 0;
xdot(61) = 0;
xdot(62) = 0;
xdot(63) = 0;
xdot(64) = 0;
xdot(65) = 0;
xdot(66) = 0;
xdot(67) = 0;
xdot(68) = 0;
xdot(69) = 0;
xdot(70) = 0;
xdot(71) = 0;
xdot(72) = 0;
xdot(73) = 0;
xdot(74) = 0;
xdot(75) = 0;
xdot(76) = 0;
xdot(77) = 0;
xdot(78) = 0;
xdot(79) = 0;


