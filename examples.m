%%  EXAMPLES
%  
%  Classical optimization problems & values etc.
%
%  BIBLIOGRAPHY:
%
%  - DEB, Kalyanmoy. "Multi-Objective optimization using evolutionary
%    algorithms". John Wiley & Sons, LTD. Kanpur, India. 2004.
%
%  - BACK, Thomas. "Evolutionary algorithms in theory and practice". Oxford
%    University Press. New York. 1996.
%
%  - BINH, Thanh; KORN, Ulrich. "MOBES: A multiobjective evolution strategy for
%    constrained optimization problems". Third International Conference on
%    Genetic Algorithms. Pag. 176-182. Czech Republic, 1997.
%
%  - BINH, Thanh. "A multiobjective evolutionary algorithm. The study cases".
%    Technical report. Barleben, Germany. 1999.

%% Beginning:
clear, clc, close all

%% Function case
fun = 1;

switch fun
  %% Constr-Ex problem (taken from DEB, pag 290.)
  case 1
    f = @(x,u) [...                       % functions to minimize
                   x(1)
                   (1 + x(2))/x(1)
               ];

    %% Constraints
    g = @(x,u) [...                       % Constraints (if they exist)
                    x(2) + 9*x(1) - 6     %  x(2) + 9*x(1) >=  6
                   -x(2) + 9*x(1) - 1     % -x(2) + 9*x(1) >=  1
                    x(1) - 0.1            %  x(1)          >=  0.1
                    1    - x(1)           %  1             >=  x(1)
                    x(2)                  %  x(2)          >=  0
                    5    - x(2)           %  5             >=  x(2)
               ];
             
    nx     = 2;                       % 'n_x' states
    nf     = length(f(nan(nx,1)));    % length of the output vector 'f(x,y)'
    ng     = length(g(nan(nx,1)));    % length of the constraint vector 'g(x,y)'
    limits = [...                     % State variables boundaries
                 0.1  1               % 0.1 <= x(1) <= 1
                 0    5               % 0   <= x(2) <= 5
             ];
    label     = {'f_{1}','f_{2}'};    % I'm going to use it in the plots

    %% Test problem BNH (Binh, Korn) - (DEB, pag. 362)
    case 2
    f = @(x,u) [...                       % functions to minimize
                   4*x(1).^2 + 4*x(2).^2
                   (x(1) - 5).^2 + (x(2) - 5).^2
               ];

    %% Constraints
    g = @(x,u) [...                                 % Constraints (if they exist)
                    25-(x(1)-5).^2 - x(2).^2        % (x(1) - 5)^2 + x(2)^2   <= 25
                    (x(1)-8).^2 + (x(2)+3).^2 - 7.7 % (x(1)-8)^2 + (x(2)+3)^2 >= 7.7
                    x(1)                            %  x(1)                   >=  0
                    5    - x(1)                     %  5                      >=  x(1)
                    x(2)                            %  x(2)                   >=  0
                    3    - x(2)                     %  3                      >=  x(2)
               ];
             
    nx     = 2;                       % 'n_x' states
    nf     = length(f(nan(nx,1)));    % length of the output vector 'f(x,y)'
    ng     = length(g(nan(nx,1)));    % length of the constraint vector 'g(x,y)'
    limits = [...                     % State variables boundaries
                 0    5               %  0 <= x(1) <= 5
                 0    3               %  0 <= x(2) <= 3
             ];
    label     = {'f_{1}','f_{2}'};    % I'm going to use it in the plots

    %% Chankong and Haimes (SRN) - (Deb, pag. 364)
    case 3
    f = @(x,u) [...                       % functions to minimize
                   (x(1) - 2).^2 + (x(2) - 1).^2 + 2
                   9*x(1) - (x(2) - 1).^2
               ];

    %% Constraints
    g = @(x,u) [...                             % Constraints (if they exist)
                    225 - x(1).^2 - x(2).^2     %  - x(1)^2 - x(2)^2 >= -225
                    -x(1) + 3*x(2) - 10         %  -x(1) + 3*x(2)    >=  10
                    x(1) + 20                   %  x(1)              >= -20
                    20   - x(1)                 %  20                >=  x(1)
                    x(2) + 20                   %  x(2)              >= -20
                    20   - x(2)                 %  20                >=  x(2)
               ];
             
    nx     = 2;                       % 'n_x' states
    nf     = length(f(nan(nx,1)));    % length of the output vector 'f(x,y)'
    ng     = length(g(nan(nx,1)));    % length of the constraint vector 'g(x,y)'
    limits = [...                     % State variables boundaries
                -20   20              % -20 <= x(1) <= 20
                -20   20              % -20 <= x(2) <= 20
             ];
    label     = {'f_{1}','f_{2}'};      % I'm going to use it in the plots

    %% Fonseca and Fleming (FON) - (DEB, pag. 339)
    case 4
    f = @(x,u) [...                       % functions to minimize
                   1 - exp(-sum((x - 1/sqrt(3)).^2))
                   1 - exp(-sum((x + 1/sqrt(3)).^2))
               ];

    %% Constraints
    g = @(x,u) [...                         % Constraints (if they exist)
                    x(1) + 4                %  x(1) >= -4
                    1    - x(1)             %  4    >=  x(1)
                    x(2) + 4                %  x(2) >= -4
                    4    - x(2)             %  4    >=  x(2)
                    x(3) + 4                %  x(3) >= -4
                    4    - x(3)             %  4    >=  x(3)
               ];
             
    nx     = 3;                       % 'n_x' states
    nf     = length(f(nan(nx,1)));    % length of the output vector 'f(x,y)'
    ng     = length(g(nan(nx,1)));    % length of the constraint vector 'g(x,y)'
    limits = [...                     % State variables boundaries
                -4    4               % -4 <= x(1) <= 4
                -4    4               % -4 <= x(2) <= 4
                -4    4               % -4 <= x(3) <= 4
             ];
    label     = {'f_{1}','f_{2}'};    % I'm going to use it in the plots

    %% Test case 4 (taken from BINH(2))
    case 5
    f = @(x,u) [...                       % functions to minimize
                   x(1).^2 - x(2)
                   -0.5*x(1) - x(2) - 1
               ];

    %% Constraints
    g = @(x,u) [...                              % Constraints (if they exist)
                    6.5 - x(1)/6 - x(2)          %  6.5 - x(1)/6 - x(2)   >= 0
                    7.5 - 0.5*x(1) - x(2)        %  7.5 - 0.5*x(1) - x(2) >= 0
                    30 - 5*x(1) - x(2)           %  30 - 5*x(1) - x(2)    >= 0
                    x(1)                         %  x(1)                  >= 0
                    x(2)                         %  x(2)                  >= 0
               ];
             
    nx     = 2;                       % 'n_x' states
    nf     = length(f(nan(nx,1)));    % length of the output vector 'f(x,y)'
    ng     = length(g(nan(nx,1)));    % length of the constraint vector 'g(x,y)'
    limits = [...                     % State variables boundaries
                -7   4                % -7 <= x(1) <= 4
                -7   4                % -7 <= x(2) <= 4
             ];
    label     = {'f_{1}','f_{2}'};       % I'm going to use it in the plots

    %% Kursawe's two-objective optimization problem (DEB, pag. 341)
    case 6
    f = @(x,u) [...                       % functions to minimize
                    -10*(exp(-0.2*sqrt(x(1)^2 + x(2)^2)) + exp(-0.2*sqrt(x(2)^2 + x(3)^2)))
                    abs(x(1)).^0.8 + abs(x(2)).^0.8 + abs(x(3)).^0.8 + 5*(sin(x(1)^3) + sin(x(2)^3) + sin(x(3)^3))
               ];

    %% Constraints
    g = @(x,u) [...                              % Constraints (if they exist)
                    x(1) + 5                     %  x(1) >= -5
                    5    - x(1)                  %  5    >= x(1)
                    x(2) + 5                     %  x(2) >= -5
                    5    - x(2)                  %  5    >= x(2)
                    x(3) + 5                     %  x(3) >= -5
                    5    - x(3)                  %  5    >= x(3)
               ];
             
    nx     = 3;                       % 'n_x' states
    nf     = length(f(nan(nx,1)));    % length of the output vector 'f(x,y)'
    ng     = length(g(nan(nx,1)));    % length of the constraint vector 'g(x,y)'
    limits = [...                     % State variables boundaries
                -5   5                % -5 <= x(1) <= 5
                -5   5                % -5 <= x(2) <= 5
                -5   5                % -5 <= x(3) <= 5
             ];
    label     = {'f_{1}','f_{2}'};      % I'm going to use it in the plots

    %% Schaffer 1 (SCH1) (DEB, pag. 338)
    case 7
    f = @(x,u) [...                       % functions to minimize
                   x.^2
                   (x-2).^2
               ];

    %% Constraints
    g = @(x,u) [...                   % Constraints (if they exist)
                    x + 5             %  x(1) >= -5
                    5 - x             %  5    >= x(1)
               ];
             
    nx     = 1;                       % 'n_x' states
    nf     = length(f(nan(nx,1)));    % length of the output vector 'f(x,y)'
    ng     = length(g(nan(nx,1)));    % length of the constraint vector 'g(x,y)'
    limits = [...                     % State variables boundaries
                -5   5                % -5 <= x(1) <= 5
             ];
    label     = {'f_{1}','f_{2}'};    % I'm going to use it in the plots

    %% (Schaffer 2 (SCH2) - (DEB, pag. 339)
    case 8
    f = @(x,u) [...                       % functions to minimize
                   -x    .*(x<=1)       + ...
                    (x-2).*(1<x & x<=3) + ...
                    (4-x).*(3<x & x<=4) + ...
                    (x-4).*(x>4);
                    (x-5).^2
               ];

    %% Constraints
    g = @(x,u) [...                   % Constraints (if they exist)
                    x  + 5            %  x(1) >= -5
                    10 - x            %  5    >= x(1)
               ];
             
    nx     = 1;                       % 'n_x' states
    nf     = length(f(nan(nx,1)));    % length of the output vector 'f(x,y)'
    ng     = length(g(nan(nx,1)));    % length of the constraint vector 'g(x,y)'
    limits = [...                     % State variables boundaries
                -5   10               % -5 <= x(1) <= 10
             ];
    label     = {'f_{1}','f_{2}'};    % I'm going to use it in the plots

    %% Poloni's two-objective problem (DEB, pag. 343)
    case 9
    A1 = 0.5*sin(1) - 2*cos(1) + sin(2)   - 1.5*cos(2);
    A2 = 1.5*sin(1) - cos(1)   + 2*sin(2) - 0.5*cos(2);
    f  = @(x,u) [...                       % functions to minimize
                    1 + (A1 - (0.5*sin(x(1)) - 2*cos(x(1)) + sin(x(2)) - 1.5*cos(x(2)))).^2 + (A2 - (1.5*sin(x(1)) - cos(x(1))   + 2*sin(x(2)) - 0.5*cos(x(2)))).^2
                    (x(1) + 3).^2 + (x(2) + 1).^2
               ];

    %% Constraints
    g = @(x,u) [...                   % Constraints (if they exist)
                    x(1) + pi         %  x(1) >= -pi
                    pi   - x(1)       %  pi   >= x(1)
                    x(2) + pi         %  x(2) >= -pi
                    pi   - x(2)       %  pi   >= x(2)
               ];
             
    nx     = 2;                       % 'n_x' states
    nf     = length(f(nan(nx,1)));    % length of the output vector 'f(x,y)'
    ng     = length(g(nan(nx,1)));    % length of the constraint vector 'g(x,y)'
    limits = [...                     % State variables boundaries
                -pi   pi              % -pi <= x(1) <= pi
                -pi   pi              % -pi <= x(2) <= pi
             ];
    label     = {'f_{1}','f_{2}'};    % I'm going to use it in the plots

    %% Viennet's two-objective problem (DEB, pag. 343)
    case 10
    f  = @(x,u) [...                       % functions to minimize
                    0.5*(x(1).^2 + x(2).^2) + sin(x(1).^2 + x(2).^2)
                    ((3*x(1) - 2*x(2) + 4).^2)/8 + ((x(1) - x(2) + 1).^2)/27 + 15
                    1/(x(1).^2 + x(2).^2 + 1) - 1.1*exp(-(x(1).^2 + x(2).^2))
               ];

    %% Constraints
    g = @(x,u) [...                    % Constraints (if they exist)
                    x(1) + 3           %  x(1) >= -3
                    3    - x(1)        %  3    >= x(1)
                    x(2) + pi          %  x(2) >= -3
                    3    - x(2)        %  3    >= x(2)
               ];
             
    nx     = 2;                        % 'n_x' states
    nf     = length(f(nan(nx,1)));     % length of the output vector 'f(x,y)'
    ng     = length(g(nan(nx,1)));     % length of the constraint vector 'g(x,y)'
    limits = [...                      % State variables boundaries
                -3   3                 % -3 <= x(1) <= 3
                -3   3                 % -3 <= x(2) <= 3
             ];
    label = {'f_{1}','f_{2}','f_{3}'}; % I'm going to use it in the plots

    %% Zitzler-Deb-Thiele's test problem 1 (ZDT1) - (DEB, pag. 356)
    case 11
    nx = 30;                         % 'n_x' states
    f  = @(x,u) [...                 % functions to minimize
                    x(1)
                    (1 + 9*sum(x(2:nx))/(nx-1))*(1 - sqrt(x(1)/(1 + 9*sum(x(2:nx))/(nx-1))))
               ];

    %% Constraints
    g = @(x,u) [...                  % Constraints (if they exist)
                    x
                    1 - x
               ];

    nf        = length(f(nan(nx,1)));   % length of the output vector 'f(x,y)'
    ng        = length(g(nan(nx,1)));   % length of the constraint vector 'g(x,y)'
    limits    = repmat([0 1], nx, 1);   % State variables boundaries
    label     = {'f_{1}','f_{2}'};      % I'm going to use it in the plots

    %% Zitzler-Deb-Thiele's test problem 2 (ZDT2) - (DEB, pag. 357)
    case 12
    nx = 30;                         % 'n_x' states
    f  = @(x,u) [...                 % functions to minimize
                    x(1)
                    (1 + 9*sum(x(2:nx))/(nx-1))*(1 - (x(1)/(1 + 9*sum(x(2:nx))/(nx-1))).^2)
               ];

    %% Constraints
    g = @(x,u) [...                  % Constraints (if they exist)
                    x
                    1 - x
               ];

    nf        = length(f(nan(nx,1)));   % length of the output vector 'f(x,y)'
    ng        = length(g(nan(nx,1)));   % length of the constraint vector 'g(x,y)'
    limits    = repmat([0 1], nx, 1);   % State variables boundaries
    label     = {'f_{1}','f_{2}'};      % I'm going to use it in the plots

    %% Zitzler-Deb-Thiele's test problem 3 (ZDT3) - (DEB, pag. 357)
    case 13
    nx = 30;                         % 'n_x' states
    f  = @(x,u) [...                 % functions to minimize
                    x(1)
                    (1 + 9*sum(x(2:nx))/(nx-1))*(1 - sqrt(x(1)/(1 + 9*sum(x(2:nx))/(nx-1))) - x(1)*sin(10*pi*x(1))/(1 + 9*sum(x(2:nx))/(nx-1)))
               ];

    %% Constraints
    g = @(x,u) [...                  % Constraints (if they exist)
                    x
                    1 - x
               ];

    nf        = length(f(nan(nx,1)));   % length of the output vector 'f(x,y)'
    ng        = length(g(nan(nx,1)));   % length of the constraint vector 'g(x,y)'
    limits    = repmat([0 1], nx, 1);   % State variables boundaries
    label     = {'f_{1}','f_{2}'};      % I'm going to use it in the plots

    %% Zitzler-Deb-Thiele's test problem 4 (ZDT4) - (DEB, pag. 358)
    case 14
    nx = 10;                         % 'n_x' states
    f  = @(x,u) [...                 % functions to minimize
                    x(1)
                    (1 + 10*(nx-1) + sum(x(2:nx).^2 - 10*cos(4*pi*x(2:nx))))*(1 - sqrt(x(1)/(1 + 10*(nx-1) + sum(x(2:nx).^2 - 10*cos(4*pi*x(2:nx))))))
               ];

    %% Constraints
    g = @(x,u) [...                  % Constraints (if they exist)
                    x(1)
                    1 - x(1)
                    x(2:end) + 5
                    5 - x(2:end)
               ];

    nf     = length(f(nan(nx,1)));   % length of the output vector 'f(x,y)'
    ng     = length(g(nan(nx,1)));   % length of the constraint vector 'g(x,y)'
    limits = [...
                0     1
                repmat([-5 5], nx-1, 1)
             ];                      % State variables boundaries
    label     = {'f_{1}','f_{2}'};   % I'm going to use it in the plots

    %% Test problem TNK - (Deb, pag. 366)
    case 15
    f  = @(x,u) [...                       % functions to minimize
                    x(1)
                    x(2)
               ];

    %% Constraints
    g = @(x,u) [...                         % Constraints (if they exist)
                    x(1).^2 + x(2).^2 - 1 - 0.1*cos(16*atan2(x(1),x(2)))  % x(1)^2 + x(2)^2 - 1 - 0.1*cos(16*atan(x(1)/x(2))) >= 0
                    0.5 - (x(1)-0.5).^2 - (x(2)-0.5).^2                   % (x(1)-0.5).^2 + (x(2)-0.5).^2                     <= 0.5
                    x(1)                                                  % x(1)                                              >= 0
                    pi   - x(1)                                           % pi                                                >= x(1)
                    x(2)                                                  % x(2)                                              >= 0
                    pi   - x(2)                                           % pi                                                >= x(2)
               ];
             
    nx     = 2;                        % 'n_x' states
    nf     = length(f(nan(nx,1)));     % length of the output vector 'f(x,y)'
    ng     = length(g(nan(nx,1)));     % length of the constraint vector 'g(x,y)'
    limits = [...                      % State variables boundaries
                 0   pi                % 0 <= x(1) <= pi
                 0   pi                % 0 <= x(2) <= pi
             ];
    label     = {'f_{1}','f_{2}'};     % I'm going to use it in the plots
    
    %% Example cantilever beam - (Deb, pag. 19)
    case 16
    rho       = 7800;                     % density (kg/m3)
    P         = 1;                        % Load (kN)
    E         = 207e6;                    % Elastic modulus (kPa)
    Sy        = 300e3;                    % Maximum stress allowed (kPa)
    delta_max = 0.005;                    % Maximum deflection allowed (mm)

    f  = @(x,u) [...                      % functions to minimize
                    0.25*rho*pi*x(2)*x(1)^2         % Weight (kg)
                    (64*P*x(2)^3)/(3*E*pi*x(1)^4)   % Deflection (m)
                ];

    %% Constraints
    g = @(x,u) [...                         % Constraints (if they exist)
                    Sy - (32*P*x(2))/(pi*x(1)^3)                          % Sy - (32*P*x(2))/(pi*x(1)^3)                 >= 0
                    delta_max - (64*P*x(2)^3)/(3*E*pi*x(1)^4)             % delta = (64*P*x(2)^3)/(3*E*pi*x(1)^4) = f(2) >= delta_max = 0.005
                    x(1) - 0.01                                           % x(1)                                         >= 10 mm
                    0.05 - x(1)                                           % 50 mm                                        >= x(1)
                    x(2) - 0.2                                            % x(2)                                         >= 200 mm
                    1    - x(2)                                           % 1000 mm                                      >= x(2)
               ];
             
    nx     = 2;                        % 'n_x' states
    nf     = length(f(nan(nx,1)));     % length of the output vector 'f(x,y)'
    ng     = length(g(nan(nx,1)));     % length of the constraint vector 'g(x,y)'
    limits = [...                      % State variables boundaries
                 0.01   0.05           % 10 mm  <= x(1) <= 50 mm
                 0.20   1.00           % 200 mm <= x(2) <= 1000 mm
             ];
    label     = {'Weight (kg)','Deflection (m)'}; % I'm going to use it in the plots

  otherwise
    error('Not supported equation');
end

%% Setting initial parameters
mu      = 100;               % parent population size
lambda  = 100;               % offspring population size
gen     = 100;               % number of generations
rec_obj = 2;                 % Type of recombination to use on object
                             % variables (Pag. 74 in (BACK))
                             % See 'recombination.m'
rec_str = 4;                 % Type of recombination to use on strategy
                             % parameters (Pag. 74 in (BACK))
u       = 0;                 % external excitation

%% Run "Elitist Non-Dominated Sorting Evolutionary Strategy" (ENSES):
[min_x, min_f] = ENSES(f, mu, lambda, gen, rec_obj, rec_str, u, nf, nx, limits, g, ng);

%% Plot initial population
for i = 1:nf-1
  for k = i+1:nf
    figure
    plot(min_f{1}(i,:),min_f{1}(k,:),'o')
    grid on
    title('Initial population, gen = 0','FontSize',16)
    xlabel(label{i},'FontSize',16);
    ylabel(label{k},'FontSize',16);
    legend('Initial population at gen = 0','location','Best')
  end
end

%% Plot function
for i = 1:nf-1
  for k = i+1:nf
    figure
    for j = 1:size(min_f,2)
      h = plot(min_f{j}(i,:),min_f{j}(k,:),'o');
      axis([min(min_f{j}(i,:)) max(min_f{j}(i,:)) min(min_f{j}(k,:)) max(min_f{j}(k,:))]);
      grid on
      title('Pareto optimal front','FontSize',18)
      xlabel(label{i},'FontSize',16);
      ylabel(label{k},'FontSize',16);
      pause(0.1)
      hold on
      if j ~= size(min_f,2)
        delete(h);
      end
    end
  end
end

if nf == 3
  figure
  plot3(min_f{end}(1,:),min_f{end}(2,:),min_f{end}(3,:),'bo')
  grid on
  title('Pareto optimal front','FontSize',18)
  xlabel(label{1},'FontSize',16);
  ylabel(label{2},'FontSize',16);
  zlabel(label{3},'FontSize',16);
end

%% END
