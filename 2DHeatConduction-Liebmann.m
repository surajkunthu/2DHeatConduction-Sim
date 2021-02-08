% Suraj Kunthu                     %%% 
% Steady State 2D Heat Conduction  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Liebmann Method
% Over-relaxation factor of Lambda, l = 1.5

% Clear all----------------------------------------------------------------
close all
clear all
clc

% Boundary Conditions------------------------------------------------------
a = 12;                     % Domain base [centimeters]
b = 6;                      % Domain height [centimeters]
%imax = 6;                  % Nodes in x-direction
imax = input ('Enter the amount of nodes in the x-direction '  );
%jmax = 3;                  % Nodes in y-direction
jmax = input ('Enter the amount of nodes in the y-direction '  );
dx = a/imax;                % x-direction step size [cm]
dy = b/jmax;                % y-direction step size [cm]
l = 1.5;                    % Over-relaxation factor
T_bottom = 0;               % [deg C]
T_left = 0;                 % [deg C]
T_right = 0;                % [deg C]
T_top = sin(0:dx:a);        % [deg C]

x = linspace(0, a, dx);     % x-vector from (0)-(a) with (dx) divisions
y = linspace(0, b, dy);     % y-vector from (0)-(b) with (dy) divisions

T = zeros(jmax+1, imax+1);  % Generate "intial guess" matrix T of i*j dimensions

T(jmax+1, 1:imax+1) = T_top;    % Top Temp. Boundary Condition
T(1, 1:imax+1) = 0;             % Bottom Temp Boundary Condition
T(1:jmax+1, 1) = 0;             % Left Temp. Boundary Condition
T(1:jmax+1, imax+1) = 0;        % Right Temp. Boundary Condition

k=0;                              % Initial Counter
%iterations = 3                   % Number of iterations
iterations = input('Enter the number of iterations '   );
for k = 1:iterations;             % Counter
   k = k+1;                       % Counter step
   To = T;                        % To = Old stored temperature per iteration
    for i = 2:(jmax+1)-1;         % Interior Nodes in x-direction
      for j = 2:(imax+1)-1;       % Interior Nodes in y-direction
          
          % Liebmann Method Finite Difference Equation
         T1(i,j) = 0.25*(T(i+1,j) + T(i-1,j) + T(i,j-1) + T(i,j+1));
          % Over-relaxation formula
         T(i,j)=l*T1(i,j)+(1-l)*T(i,j);
      end
    end
end

T(3, 1:(imax+1))
Total_Nodes_in_Mesh = (imax+1)*(jmax+1)

% Graph Setup
sizex=0:dx:a;
sizey=0:dy:b;
figure(1)
contourf(sizex,sizey,T)
colorbar
xlabel({'X [cm]'});
ylabel({'Y [cm]'});
title('Over-relaxed Steady State 2D Heat Conduction Temperature Distribution')
annotation('textbox',[0.85, 0.95, 0.1, 0.03],'String',{'Temp. ^{\circ}C'});

figure(2)
contour3(sizex,sizey,T); hold on; surf(sizex,sizey,T, 'Edgecolor', 'none');
colorbar
xlabel({'X [cm]'});
ylabel({'Y [cm]'});
zlabel({'Temperature [deg C]'})
title('Over-relaxed Steady State 2D Heat Conduction Temperature Distribution')
annotation('textbox',[0.85, 0.95, 0.1, 0.03],'String',{'Temp. ^{\circ}C'});

figure(3)
plot(sizex, T)%T(3, 1:(imax+1)))
xlabel({'X [cm]'});
ylabel({'T [deg C]'});
title('Over-relaxed Steady State 2D Heat Conduction Temperature Profiles')
legend('Nodal Temperatures [deg C]')
