

figure
subplot(4,2,1);
plot (x_v(:,1),'b')
legend('Glucose 6P in fuzzy model');
ylabel('Activity/Concentration')
xlabel('Time')

subplot(4,2,2);
plot (hypoxia(:,2),'m')
legend('Glucose 6P in Hypoxia result');
ylabel('Activity/Concentration')
xlabel('Time')


subplot(4,2,3);
plot (x_v(:,3),'b')
legend('Fructose 6P in fuzzy model');
ylabel('Activity/Concentration')
xlabel('Time')

subplot(4,2,4);
plot (hypoxia(:,4),'m')
legend('Fructose 6P in Hypoxia result');
ylabel('Activity/Concentration')
xlabel('Time')


subplot(4,2,5);
plot (x_v(:,4),'b')
legend('Fructose 1,6BP in fuzzy model');
ylabel('Activity/Concentration')
xlabel('Time')

subplot(4,2,6);
plot (hypoxia(:,5),'m')
legend('Fructose 1,6BP in Hypoxia result');
ylabel('Activity/Concentration')
xlabel('Time')


subplot(4,2,7);
plot (x_v(:,6),'b')
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
plot (x_v(:,9),'b')
legend('3PG in fuzzy model');
ylabel('Activity/Concentration')
xlabel('Time')

subplot(4,2,2);
plot (hypoxia(:,10),'m')
legend('3PG in Hypoxia result');
ylabel('Activity/Concentration')
xlabel('Time')


subplot(4,2,3);
plot (x_v(:,11),'b')
legend('PEP in fuzzy model');
ylabel('Activity/Concentration')
xlabel('Time')

subplot(4,2,4);
plot (hypoxia(:,12),'m')
legend('PEP in Hypoxia result');
ylabel('Activity/Concentration')
xlabel('Time')


subplot(4,2,5);
plot (x_v(:,12),'b')
legend('Pyruvate in fuzzy model');
ylabel('Activity/Concentration')
xlabel('Time')

subplot(4,2,6);
plot (hypoxia(:,13),'m')
legend('Pyruvate in Hypoxia result');
ylabel('Activity/Concentration')
xlabel('Time')


subplot(4,2,7);
plot (x_v(:,13),'b')
legend('LACTATE in fuzzy model');
ylabel('Activity/Concentration')
xlabel('Time')

subplot(4,2,8);
plot (hypoxia(:,14),'m')
legend('LACTATE in Hypoxia result');
ylabel('Activity/Concentration')
xlabel('Time')




