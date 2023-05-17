function PlotAircraftSim(TOUT, aircraft_state, control_surfaces, wind_inertial, maintitle,i,col)



len_sim = length(TOUT);
    
    
    for j = 1:len_sim
        [flight_angles] = FlightPathAnglesFromState(aircraft_state(j,:)');
        chi(j,1) = (180/pi)*flight_angles(2,1);
        gamma(j,1) = (180/pi)*flight_angles(3,1);
        
        wind_body = TransformFromInertialToBody(wind_inertial, aircraft_state(j,4:6)');
        air_rel_vel_body = aircraft_state(j,7:9)' - wind_body;

        wind_angles(:,j) = WindAnglesFromVelocityBody(air_rel_vel_body);
    end


    
    

%%%%%%%%%%%%%%%%%%%%%%%%
figure(8*i - 7);
sgtitle(maintitle);
subplot(311);
h1= plot(TOUT, aircraft_state(:,1),col);hold on;
title('Position v Time');   
ylabel('X [m]')   
grid on
legend('No Control','Control')

subplot(312);
plot(TOUT, aircraft_state(:,2),col);hold on;
 ylabel('Y [m]')
 grid on
 legend('No Control','Control')
 
subplot(313);
plot(TOUT, aircraft_state(:,3),col);hold on;
ylabel('Z [m]')    
xlabel('time [sec]');
grid on
legend('No Control','Control')

%%%%%%%%%%%%%%%%%%%%%%%%
figure(8*i - 6);
sgtitle(maintitle);
subplot(311);
plot(TOUT, (180/pi)*aircraft_state(:,4),col);hold on;
title('Euler Angles v Time');   
ylabel('Roll [deg]')    
grid on
legend('No Control','Control')

subplot(312);
plot(TOUT, (180/pi)*aircraft_state(:,5),col);hold on;
 ylabel('Pitch [deg]')   
 grid on
legend('No Control','Control')
 
subplot(313);
plot(TOUT, (180/pi)*aircraft_state(:,6),col);hold on;
ylabel('Yaw [deg]')    
xlabel('time [sec]');
grid on
legend('No Control','Control')

%%%%%%%%%%%%%%%%%%%%%%%%
figure(8*i - 5);
sgtitle(maintitle);
subplot(311);
plot(TOUT, aircraft_state(:,7),col);hold on;
title('Velocity v Time');   
ylabel('uE [m/s]')    
grid on
legend('No Control','Control')

subplot(312);
plot(TOUT, aircraft_state(:,8),col);hold on;
 ylabel('vE [m/s]')    
 grid on
legend('No Control','Control')
 
subplot(313);
plot(TOUT, aircraft_state(:,9),col);hold on;
ylabel('wE [m/s]')    
xlabel('time [sec]');
grid on
legend('No Control','Control')

%%%%%%%%%%%%%%%%%%%%%%%%
figure(8*i - 4);
sgtitle(maintitle);
subplot(311);
plot(TOUT, (180/pi)*aircraft_state(:,10),col);hold on;
title('Angular Velocity v Time');   
ylabel('p [deg/s]')    
grid on
legend('No Control','Control')

subplot(312);
plot(TOUT, (180/pi)*aircraft_state(:,11),col);hold on;
 ylabel('q [deg/s]') 
 grid on
legend('No Control','Control')
 
subplot(313);
plot(TOUT, (180/pi)*aircraft_state(:,12),col);hold on;
ylabel('r [deg/s]')    
xlabel('time [sec]');
grid on
legend('No Control','Control')


%%%%%%%%%%%%%%%%%%%%%%%%
figure(8*i - 3);
sgtitle(maintitle);
plot3(aircraft_state(:,1),aircraft_state(:,2),-aircraft_state(:,3),col);hold on;
grid on
legend('No Control','Control')

%%%%%%%%%%%%%%%%%%%%%%%%
if (~isempty(control_surfaces))
    figure(8*i - 2);
    sgtitle(maintitle);
    subplot(411);
    plot(TOUT, control_surfaces(:,1)*(180/pi),col);hold on;
    title('Control Surfaces v Time');   
    ylabel('\delta_e [deg]')    
    grid on
    legend('No Control','Control')
    
    subplot(412);
    plot(TOUT, control_surfaces(:,2)*(180/pi),col);hold on;
    ylabel('\delta_a [deg]')      
    grid on
    legend('No Control','Control')
    
    subplot(413);
    plot(TOUT, control_surfaces(:,3)*(180/pi),col);hold on;
    ylabel('\delta_r [deg]')   
    grid on
    legend('No Control','Control')
    
    subplot(414);
    plot(TOUT, control_surfaces(:,4),col);hold on;
    ylabel('\delta_t [frac]')     
    xlabel('time [sec]');
    grid on
    legend('No Control','Control')

end

%%%%%%%%%%%%%%%%%%%%%%%%
figure(8*i - 1);
sgtitle(maintitle);
subplot(311);
h7= plot(TOUT, wind_angles(1,:) ,col);hold on;
title('Wind Angles v Time');   
ylabel('V_a [m/s]')   
grid on
legend('No Control','Control')

subplot(312);
plot(TOUT, (180/pi)*wind_angles(2,:),col);hold on;
ylabel('\beta [deg]')  
grid on
legend('No Control','Control')
 
subplot(313);
plot(TOUT, (180/pi)*wind_angles(3,:),col);hold on;
ylabel('\alpha [deg]')    
xlabel('time [sec]');
grid on
legend('No Control','Control')

figure(8*i);
sgtitle(maintitle);
subplot(211)
plot(TOUT, chi, col);hold on;
ylabel('Course angle \chi [deg]')   
grid on
legend('No Control','Control')

subplot(212)
plot(TOUT, gamma, col);hold on;
ylabel('Flight path angle \gamma [deg]') 
grid on
legend('No Control','Control')


