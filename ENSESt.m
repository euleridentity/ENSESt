function [min_x, min_f] = ENSESt(f, mu, lambda, gen, rec_obj, rec_str, u, nf, n, limits, g, ng, x_0, sigma)
%%  [min_x, min_f] = ENSESt(f, mu, lambda, gen, rec_obj, rec_str, u, nf, n, limits, g, ng, x_0, sigma)
%
%  This function implements the "Elitist Non-Dominated Sorting Evolution
%  Strategy" (ENSES) which is an elitist multi-objective evolutionary
%  algorithm. The algorithm is based in NSGA-II, which can be studied in
%  section 6.2 of (DEB), but it uses "Evolution strategy" as
%  evolutionary algorithm.
%
%  INPUT DATA:
%
%  - f:       Objective function (handle function: f(x,u))
%  - mu:      Parent population size (positive integer number)
%  - lambda:  Offspring population size (positive integer number)
%  - gen:     Number of generations (positive integer number)
%  - rec_obj: Type of recombination to use on objective variables
%             (Pag. 74 in (BACK)):
%                  * 1   = No recombination
%                  * 2   = Discrete recombination
%                  * 3   = Panmictic discrete recombination
%                  * 4   = Intermediate recombination
%                  * 5   = Panmictic intermediate recombination
%                  * 6   = Generalized intermediate recombination
%                  * 7   = Panmictic generalized intermediate recombination
%  - rec_str: Type of recombination to use on strategy parameters
%             (Pag. 74 in (BACK)).
%  - u:       External excitation (if it does not exist, type 0 (zero))
%  - nf:      Length of the handle function vector (length(f) x 1 vector)
%  - n:       Length of the vector x_0 (positive integer number)
%  - limits:  Matrix with the limits of the variables (nx2 matrix). The
%             first column is the lower boundary, the second column is the
%             upper boundary.
%  - g:       Constraints (handle function g(x,u))
%  - ng:      Length of the constraint handle function vector
%             (length(g) x 1 vector)
%  - x_0:     Starting point (optional) (nxmu matrix)
%  - sigma:   Cell with covariance matrices (Optional) (1 x mu cell; each
%             cell has to have an nxn symmetric matrix)
%
%  OUTPUT DATA:
%
%  - min_x: Vector which minimizes the objective function 'f' (nx1 vector)
%  - min_f: Achived minimum of the function 'f' (length(f) x 1 vector)
%
%  NOTES:
%
%  - The implementation of this multi-objective optimization technique uses
%    "evolution strategies" as evolutionary algorithm.
%  - The algorithm only works with the 'mu + lambda' selection scheme.
%  - It only works in minimization of the objective functions. Because of
%    the 'Duality principle', you can transform the original objective for
%    maximization into an objective for minimization. This can be achieved
%    just multiplying the objective function by -1.
%
%  BIBLIOGRAPHY:
%
%  - DEB, Kalyanmoy. "Multi-Objective optimization using evolutionary
%    algorithms". John Wiley & Sons, LTD. Kanpur, India. 2004.
%
%  - BACK, Thomas. "Evolutionary algorithms in theory and practice". Oxford
%    University Press. New York. 1996.

%% Beggining
ENSESt('Examples', 100, 25, 1000, 2, 1, 0, length(f)*1, 10, [-100,100], g, ng, x_0, sigma)

f='Examples';
mu=100;
lamda=25;
gen=1000;
rec_obj=2;
rec_str=1;



%% Initialization:
if exist('x_0','var')
  if exist('sigma','var')
    [x_0, sigma, alpha] = validate_sizes_es(mu, n, limits, x_0, sigma);
  else
    [x_0, sigma, alpha] = validate_sizes_es(mu, n, limits, x_0);
  end
else
  [x_0, sigma, alpha] = validate_sizes_es(mu, n, limits);
end

min_x = cell(1,gen);                  % allocate space in memory for min_x
min_f = cell(1,gen);                  % allocate space in memory for min_f

min_x{1} = x_0;                       % first point
value    = zeros(nf,mu);              % allocate space for function evaluation
for i = 1:mu
  value(:,i) = f(x_0(:,i),u);
end
min_f{1} = value;                     % first approximation
j        = 1;                         % generations counter

%% Start ENSESt
while (j < gen)
  %% Print report:
  msg = fprintf('\tGeneration j = %4d\n',j);
  
  %% Recombine:
  [xr,sigmar,alphar] = recombination(rec_obj,rec_str,n,mu,lambda,min_x{j},sigma,alpha);

  %% Mutation:
  [xm,sigmam,alpham] = mutation_es(n,lambda,xr,sigmar,alphar,limits);
  
  %% Elitism
  P = [min_x{j} xm];                  % 'mu + lambda'
  
  %% Evaluation:
  ml   = mu + lambda;                 % mu + lambda = size(P,2)
  phie = zeros(nf,ml);
  for i = 1:ml
    phie(:,i) = f(P(:,i),u);
  end
  
  %% Identify fronts and carry out crowded tournament selection:
  [min_x{j+1},sigma,alpha] = non_dominated_sorting_algorithm_es(phie,P,u,sigmam,sigma,alpham,alpha,mu,lambda,nf,n,g,ng);

  %% Store better results:
  value = zeros(nf,mu);               % allocate space for function evaluation
  for i = 1:mu
    value(:,i) = f(min_x{j+1}(:,i),u);
  end
  min_f{j+1} = value;                 % next approximation
  
  %% Increase generation counter:
  j = j+1;

  fprintf(repmat('\b',1,msg));
end

end
%% END
