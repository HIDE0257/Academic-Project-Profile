function PlotAircraftSim_MyVersion(TOUT, aircraft_state, control_surfaces, background_wind_array, col)


%%%%%%%%%%%%%%%%%%%%%%%%
figure(1);
sgtitle(sprintf('Position vs. Time with a Wind of [%d; %d; %d] m/s', background_wind_array))
subplot(311);
h1= plot(TOUT, aircraft_state(:,1),col,'LineWidth',1);hold on;   
ylabel('X [m]')    
xlabel('time [sec]')
legend('Non-Linear','Linearized')

subplot(312);
plot(TOUT, aircraft_state(:,2),col,'LineWidth',1);hold on;
ylabel('Y [m]')    
xlabel('time [sec]')
legend('Non-Linear','Linearized')

subplot(313);
plot(TOUT, aircraft_state(:,3),col,'LineWidth',1);hold on;
ylabel('Z [m]')    
xlabel('time [sec]')
legend('Non-Linear','Linearized')

%%%%%%%%%%%%%%%%%%%%%%%%
figure(2);
sgtitle(sprintf('Eular Angles vs. Time with a Wind of [%d; %d; %d] m/s', background_wind_array))
subplot(311);
plot(TOUT, (180/pi)*aircraft_state(:,4),col,'LineWidth',1);hold on;
ylabel('Roll \phi [deg]')    
xlabel('Time [sec]')
legend('Non-Linear','Linearized')

subplot(312);
plot(TOUT, (180/pi)*aircraft_state(:,5),col,'LineWidth',1);hold on;
ylabel('Pitch \theta [deg]')   
xlabel('Time [sec]')
legend('Non-Linear','Linearized')

subplot(313);
plot(TOUT, (180/pi)*aircraft_state(:,6),col,'LineWidth',1);hold on;
ylabel('Yaw \psi [deg]')    
xlabel('time [sec]');
legend('Non-Linear','Linearized')

%%%%%%%%%%%%%%%%%%%%%%%%
figure(3);
sgtitle(sprintf('Inertial Velovity vs. Time with a Wind of [%d; %d; %d] m/s', background_wind_array))

subplot(311);
plot(TOUT, aircraft_state(:,7),col,'LineWidth',1);hold on;
ylabel('u^{E} [m/s]')    
xlabel('Time [sec]');
legend('Non-Linear','Linearized')

subplot(312);
plot(TOUT, aircraft_state(:,8),col,'LineWidth',1);hold on;
 ylabel('v^{E} [m/s]')    
 xlabel('Time [sec]');
 legend('Non-Linear','Linearized')

subplot(313);
plot(TOUT, aircraft_state(:,9),col,'LineWidth',1);hold on;
ylabel('w^{E} [m/s]')    
xlabel('Time [sec]');
legend('Non-Linear','Linearized')

%%%%%%%%%%%%%%%%%%%%%%%%
figure(4);
sgtitle(sprintf('Angular Velovity vs. Time with a Wind of [%d; %d; %d] m/s', background_wind_array))

subplot(311);
plot(TOUT, (180/pi)*aircraft_state(:,10),col,'LineWidth',1);hold on;  
ylabel('p [deg/s]')    
xlabel('Time [sec]');
legend('Non-Linear','Linearized')

subplot(312);
plot(TOUT, (180/pi)*aircraft_state(:,11),col,'LineWidth',1);hold on;
 ylabel('q [deg/s]')   
 xlabel('Time [sec]');
 legend('Non-Linear','Linearized')

subplot(313);
plot(TOUT, (180/pi)*aircraft_state(:,12),col,'LineWidth',1);hold on;
ylabel('r [deg/s]')    
xlabel('Time [sec]');
legend('Non-Linear','Linearized')

%%%%%%%%%%%%%%%%%%%%%%%%
figure(5);
title(sprintf('Trajectory over Time with a Wind of [%d; %d; %d] m/s', background_wind_array))
plot3(aircraft_state(:,1),aircraft_state(:,2),-aircraft_state(:,3),col,'LineWidth',2);hold on
xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')
grid on
legend('Non-Linear','Linearized')
% axis([0 600 100 500 2436 2440])

%%%%%%%%%%%%%%%%%%%%%%%%
if (~isempty(control_surfaces))
   figure(6);
  sgtitle(sprintf('Control Surface vs. Time with a Wind of [%d; %d; %d] m/s', background_wind_array))

   subplot(411);
   plot(TOUT, control_surfaces(:,1),col,'LineWidth',1);hold on;
   ylabel('Elevator [rad]')  
   xlabel('Time [sec]');
   legend('Non-Linear','Linearized')

   subplot(412);
   plot(TOUT, control_surfaces(:,2),col,'LineWidth',1);hold on;
   ylabel('Aileron [rad]') 
   xlabel('Time [sec]');
   legend('Non-Linear','Linearized')

   subplot(413);
   plot(TOUT, control_surfaces(:,3),col,'LineWidth',1);hold on;
   ylabel('Rudder [rad]')      
  xlabel('Time [sec]');
  legend('Non-Linear','Linearized')

  
   subplot(414);
   plot(TOUT, control_surfaces(:,4),col,'LineWidth',1);hold on;
   ylabel('Throttle [frac]')    
   xlabel('Time [sec]');
   legend('Non-Linear','Linearized')

end


euler_angles = aircraft_state(:,4:6);
vel_body = aircraft_state(:,7:9);

for i = 1:length(vel_body)
    wind_body(:,i) = TransformFromInertialToBody(background_wind_array, euler_angles(i,:)');
    vel_air_relative(:,i) = vel_body(i,:)' - wind_body(:,i);
    Va(i) = norm(vel_air_relative(:,i));
    beta(i) = asin(vel_body(i,2)/Va(i));
    alpha(i) = atan(vel_body(i,3)/vel_body(i,1));
end

figure(7)
sgtitle(sprintf('Airspeed and Wind Angles over Time with a Wind of [%d; %d; %d] m/s', background_wind_array))

subplot(311);
plot(TOUT, Va,col,'LineWidth',1);hold on;
ylabel('Air Speed [m/s]')   
xlabel('Time [sec]');
legend('Non-Linear','Linearized')


subplot(312);
plot(TOUT, beta,col,'LineWidth',1);hold on;
ylabel('Side Slip Angle [rad]') 
xlabel('Time [sec]');
legend('Non-Linear','Linearized')


subplot(313);
plot(TOUT, alpha,col,'LineWidth',1);hold on;
ylabel('Angle of Atack [rad]')   
xlabel('Time [sec]');
legend('Non-Linear','Linearized')

