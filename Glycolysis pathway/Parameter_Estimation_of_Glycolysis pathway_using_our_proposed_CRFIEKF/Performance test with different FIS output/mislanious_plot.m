
% load('ER_stress_Insulin_Cardio__Pathway_Normal_condition.mat');
% load('ER_stress_Insulin_Cardio__Pathway_Stressed_condition.mat');

% figure
% % x = linspace(1,70001);
% plot (xe(3,:))
% hold on
% plot (observed(3,:))
% legend('estimated','observed');
% % axis([0 70001 .48 .56]);
% ylabel('Activity/Concentration')
% xlabel('Time')


figure
subplot(4,2,1);
plot (x_v(:,2),'b')
legend('Glucose 6P in fuzzy model');
ylabel('Activity/Concentration')
xlabel('Time')

subplot(4,2,2);
plot (hypoxia(:,2),'m')
legend('Glucose 6P in Hypoxia result');
ylabel('Activity/Concentration')
xlabel('Time')


subplot(4,2,3);
plot (x_v(:,4),'b')
legend('Fructose 6P in fuzzy model');
ylabel('Activity/Concentration')
xlabel('Time')

subplot(4,2,4);
plot (hypoxia(:,4),'m')
legend('Fructose 6P in Hypoxia result');
ylabel('Activity/Concentration')
xlabel('Time')


subplot(4,2,5);
plot (x_v(:,5),'b')
legend('Fructose 1,6BP in fuzzy model');
ylabel('Activity/Concentration')
xlabel('Time')

subplot(4,2,6);
plot (hypoxia(:,5),'m')
legend('Fructose 1,6BP in Hypoxia result');
ylabel('Activity/Concentration')
xlabel('Time')


subplot(4,2,7);
plot (x_v(:,7),'b')
legend('DHAP in fuzzy model');
ylabel('Activity/Concentration')
xlabel('Time')

subplot(4,2,8);
plot (hypoxia(:,7),'m')
legend('DHAP in Hypoxia result');
ylabel('Activity/Concentration')
xlabel('Time')







figure
subplot(4,2,1);
plot (x_v(:,10),'b')
legend('3PG in fuzzy model');
ylabel('Activity/Concentration')
xlabel('Time')

subplot(4,2,2);
plot (hypoxia(:,10),'m')
legend('3PG in Hypoxia result');
ylabel('Activity/Concentration')
xlabel('Time')


subplot(4,2,3);
plot (x_v(:,12),'b')
legend('PEP in fuzzy model');
ylabel('Activity/Concentration')
xlabel('Time')

subplot(4,2,4);
plot (hypoxia(:,12),'m')
legend('PEP in Hypoxia result');
ylabel('Activity/Concentration')
xlabel('Time')


subplot(4,2,5);
plot (x_v(:,13),'b')
legend('Pyruvate in fuzzy model');
ylabel('Activity/Concentration')
xlabel('Time')

subplot(4,2,6);
plot (hypoxia(:,13),'m')
legend('Pyruvate in Hypoxia result');
ylabel('Activity/Concentration')
xlabel('Time')


subplot(4,2,7);
plot (x_v(:,14),'b')
legend('LACTATE in fuzzy model');
ylabel('Activity/Concentration')
xlabel('Time')

subplot(4,2,8);
plot (hypoxia(:,14),'m')
legend('LACTATE in Hypoxia result');
ylabel('Activity/Concentration')
xlabel('Time')





% % % 
% % % figure
% % % subplot(3,3,1);
% % % plot (x_v(:,2))
% % % hold on
% % % plot (hypoxia(:,2))
% % % legend('Hypoxia result','Result after fuzzy PM');
% % % % axis([0 10000 .4 .6]);
% % % ylabel('Activity/Concentration')
% % % xlabel('Time')
% % % 
% % % subplot(3,3,2);
% % % plot (x_v(:,4))
% % % hold on
% % % plot (hypoxia(:,4))
% % % legend('Hypoxia result','Result after fuzzy PM');
% % % % axis([0 10000 .4 .6]);
% % % ylabel('Activity/Concentration')
% % % xlabel('Time')
% % % 
% % % subplot(3,3,3);
% % % plot (x_v(:,5))
% % % hold on
% % % plot (hypoxia(:,5))
% % % legend('Hypoxia result','Result after fuzzy PM');
% % % % axis([0 10000 .2 .8]);
% % % ylabel('Activity/Concentration')
% % % xlabel('Time')
% % % 
% % % subplot(3,3,4);
% % % plot (x_v(:,7))
% % % hold on
% % % plot (hypoxia(:,7))
% % % legend('Hypoxia result','Result after fuzzy PM');
% % % % axis([0 10000 .4 .8]);
% % % ylabel('Activity/Concentration')
% % % xlabel('Time')
% % % 
% % % subplot(3,3,5);
% % % plot (x_v(:,10))
% % % hold on
% % % plot (hypoxia(:,10))
% % % legend('Hypoxia result','Result after fuzzy PM');
% % % % axis([0 10000 .4 .8]);
% % % ylabel('Activity/Concentration')
% % % xlabel('Time')
% % % 
% % % subplot(3,3,6);
% % % plot (x_v(:,12))
% % % hold on
% % % plot (hypoxia(:,12))
% % % legend('Hypoxia result','Result after fuzzy PM');
% % % % axis([0 10000 .4 .8]);
% % % ylabel('Activity/Concentration')
% % % xlabel('Time')
% % % 
% % % 
% % % subplot(3,3,7);
% % % plot (x_v(:,13))
% % % hold on
% % % plot (hypoxia(:,13))
% % % legend('Hypoxia result','Result after fuzzy PM');
% % % % axis([0 10000 .4 .8]);
% % % ylabel('Activity/Concentration')
% % % xlabel('Time')
% % % 
% % % subplot(3,3,8);
% % % plot (x_v(:,14))
% % % hold on
% % % plot (hypoxia(:,14))
% % % legend('Hypoxia result','Result after fuzzy PM');
% % % % axis([0 10000 .4 .8]);
% % % ylabel('Activity/Concentration')
% % % xlabel('Time')
% % % 
% % % subplot(3,3,9);
% % % plot (x_v(:,9))
% % % hold on
% % % plot (hypoxia(:,9))
% % % legend('Hypoxia result','Result after fuzzy PM');
% % % % axis([0 10000 .4 .8]);
% % % ylabel('Activity/Concentration')
% % % xlabel('Time')
% % % 
% % % 
% % % 
% % % 
% % % 
% % % 
% % % 
% % % figure
% % % subplot(3,3,1);
% % % plot (x_v(:,1))
% % % hold on
% % % plot (hypoxia(:,1))
% % % legend('Hypoxia result','Result after fuzzy PM');
% % % % axis([0 10000 .4 .6]);
% % % ylabel('Activity/Concentration')
% % % xlabel('Time')
% % % 
% % % subplot(3,3,2);
% % % plot (x_v(:,3))
% % % hold on
% % % plot (hypoxia(:,3))
% % % legend('Hypoxia result','Result after fuzzy PM');
% % % % axis([0 10000 .4 .6]);
% % % ylabel('Activity/Concentration')
% % % xlabel('Time')
% % % 
% % % subplot(3,3,3);
% % % plot (x_v(:,6))
% % % hold on
% % % plot (hypoxia(:,6))
% % % legend('Hypoxia result','Result after fuzzy PM');
% % % % axis([0 10000 .2 .8]);
% % % ylabel('Activity/Concentration')
% % % xlabel('Time')
% % % 
% % % subplot(3,3,4);
% % % plot (x_v(:,8))
% % % hold on
% % % plot (hypoxia(:,8))
% % % legend('Hypoxia result','Result after fuzzy PM');
% % % % axis([0 10000 .4 .8]);
% % % ylabel('Activity/Concentration')
% % % xlabel('Time')
% % % 
% % % subplot(3,3,5);
% % % plot (x_v(:,9))
% % % hold on
% % % plot (hypoxia(:,9))
% % % legend('Hypoxia result','Result after fuzzy PM');
% % % % axis([0 10000 .4 .8]);
% % % ylabel('Activity/Concentration')
% % % xlabel('Time')
% % % 
% % % subplot(3,3,6);
% % % plot (x_v(:,15))
% % % hold on
% % % plot (hypoxia(:,15))
% % % legend('Hypoxia result','Result after fuzzy PM');
% % % % axis([0 10000 .4 .8]);
% % % ylabel('Activity/Concentration')
% % % xlabel('Time')
% % % 
% % % 
% % % subplot(3,3,7);
% % % plot (x_v(:,16))
% % % hold on
% % % plot (hypoxia(:,16))
% % % legend('Hypoxia result','Result after fuzzy PM');
% % % % axis([0 10000 .4 .8]);
% % % ylabel('Activity/Concentration')
% % % xlabel('Time')
% % % 
% % % subplot(3,3,8);
% % % plot (x_v(:,17))
% % % hold on
% % % plot (hypoxia(:,17))
% % % legend('Hypoxia result','Result after fuzzy PM');
% % % % axis([0 10000 .4 .8]);
% % % ylabel('Activity/Concentration')
% % % xlabel('Time')
% % % 
% % % subplot(3,3,9);
% % % plot (x_v(:,18))
% % % hold on
% % % plot (hypoxia(:,18))
% % % legend('Hypoxia result','Result after fuzzy PM');
% % % % axis([0 10000 .4 .8]);
% % % ylabel('Activity/Concentration')
% % % xlabel('Time')
% % % 
