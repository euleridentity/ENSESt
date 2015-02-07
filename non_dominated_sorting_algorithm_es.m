function [min_x,min_sigma,min_alpha] = non_dominated_sorting_algorithm_es(phi,P,u,sigma_o,sigma_p,alpha_o,alpha_p,mu,lambda,n_f,n_x,g,ng)
%% [min_x,min_sigma,min_alpha] = non_dominated_sorting_algorithm_es(phi,P,sigma_o,u,sigma_p,alpha_o,alpha_p,mu,lambda,n_f,n_x,g,ng)
%  This function identify the different fronts using an O(MN^2) non-dominated
%  sorting procedure. Taken from section 2.4.7 of the bibliography given
%  below. It also carries out the "Crowded Tournament Selection", so the
%  function returns the parent population of the next generation.
%
%  INPUT DATA:
%
%  - phi:     Matrix with the output of the functions evaluated at 'x'
%             (nf x (mu + lambda) matrix)
%  - P:       Matrix with the parents and offspring population
%             (n x (mu + lambda) matrix)
%  - u:      External excitation (if it does not exist, type 0 (zero))
%  - sigma_o: Offspring covariances matrices (cell 1 x lambda)
%  - sigma_p: Parents covariances matrices (cell 1 x mu)
%  - alpha_o: Offspring rotation angles of the covariance matrices 'sigma_o'
%             (cell 1 x lambda)
%  - alpha_p: Parents rotation angles of the covariance matrices 'sigma_o'
%             (cell 1 x mu)
%  - mu:      Parent population size
%  - lambda:  Offspring population size
%  - n_f:     Number of functions to optimize (Number of rows of 'phi')
%  - n_x:     Number of state variables
%  - g:       Constrain handle function (ng x 1 vector)
%  - ng:      Number of constraints (positive integer number)
%
%  OUTPUT DATA:
%
%  - min_x:     Matrix with the parents of the next generation (n_x x mu matrix)
%  - min_sigma: Covariances matrices of the parents of the next genration
%               (cell 1 x mu). Each cell has an nxn symmetric matrix
%  - min_alpha: Rotation angles of the covariances matrices of the parents of
%               the next generation (cell 1 x mu). Each cell has an 
%               (n x n) matrix
%
%  BIBLIOGRAPHY:
%
%  - DEB, Kalyanmoy. "Multi-Objective optimization using evolutionary
%    algorithms". John Wiley & Sons, LTD. Kanpur, India. 2004.


%% Beginning

%%  Population size
m = mu + lambda;

%%  Non-domination Sorting (nds) structure:
%   nds.n_i: Number of individuals which dominate individual 'i'
%   nds.S_i: Set of individuals that individual 'i' dominates
nds = repmat(struct('n_i',0, 'S_i', []),[1,m]);

%%  Auxiliary structure, where it is going to be stored the rank of each
%   member of the population and their distances in the 'Crowding distance'
%   procedure:
%   pop.rank:     Rank of each individual
pop = repmat(struct('rank',[]),[1,m]);
for i = 1:m
    pop(i).rank = 0;
end

%% Evaluation of constraints:
phi_r = zeros(ng,m);                % Vector with the constraint values of each
                                    % member of the population
for i = 1:m
  phi_r(:,i) = g(P(:,i),u);
end

%% Check feasibility of each individual:
n_viol   = zeros(1,m);
sum_viol = zeros(1,m);
for i = 1:m
  for j = 1:ng
    if phi_r(j,i) < 0
      n_viol(i)   = n_viol(i) + 1;
      sum_viol(i) = sum_viol(i) + abs(phi_r(j,i));
    end
  end
end

%%  Domination matrix for efficiency
dom_mat  = domination_matrix(n_viol, sum_viol, phi, m, n_f);

%%  Fast non-dominated sorting procedure: Compute 'n_i' and 'S_i' of each
%   indivudal
for i = 1:m-1
  for j = i+1:m
    if(dom_mat(i,j) == 1)               % 'i' dominates 'j'
      nds(j).n_i = nds(j).n_i + 1;
      nds(i).S_i = [nds(i).S_i  j];
    elseif(dom_mat(i,j) == -1)          % 'j' dominates 'i'
      nds(i).n_i = nds(i).n_i + 1;
      nds(j).S_i = [nds(j).S_i  i];
    end
  end
end

%%  Identify first front:
front(1).index = [];                    % Initialize variable
for i = 1:m
  if nds(i).n_i == 0
    pop(i).rank = 1;
    front(1).index = [front(1).index  i];
  end
end

%%  Identify subsequent fronts
nF = 1;
while ~isempty(front(nF).index)
    Q = [];
    for i = front(nF).index
      for j = nds(i).S_i
        nds(j).n_i = nds(j).n_i - 1;
          if nds(j).n_i == 0
            pop(j).rank = nF + 1;
            Q = [Q  j];
          end
      end
    end
    nF = nF + 1;
    front(nF).index = Q;
end
front(nF) = [];                     % delete the last empty front set

%% Arrange elements of the first front (output of the function)
Fr       = cell(1,nF-1);
sigma    = [sigma_p sigma_o];
alpha    = [alpha_p alpha_o];
sigma_or = cell(1,m);
alpha_or = cell(1,m);
kk = 1;
for i = 1:nF-1
  idx = front(i).index;
  Fr{i} = P(:,idx);
  for j = 1:length(idx)
    sigma_or{kk} = sigma{idx(j)};   % Re-arrange covariance matrices
    alpha_or{kk} = alpha{idx(j)};   % Re-arrange rotation angles of the
                                    % covariance matrices
    kk = kk + 1;
  end
  
end

%% Crowded tournament selection
total_size = 0;                     % Initialize size counter
nF         = 1;                     % Front counter
min_x      = zeros(n_x,mu);         % initialize 'min_x'
min_sigma  = cell(1,mu);            % initialize 'min_sigma'
min_alpha  = cell(1,mu);            % initialize 'min_alpha'
kk         = 1;                     % auxiliary variable (counter)
while (total_size <= mu)
  tmp = size(Fr{nF},2);
  total_size = total_size + tmp;
  if total_size <= mu
    for k = 1:size(Fr{nF},2)
      min_x(:,kk)   = Fr{nF}(:,k);
      min_sigma{kk} = sigma_or{kk};
      min_alpha{kk} = alpha_or{kk};
      kk            = kk + 1;
    end
    nF = nF + 1;                    % Increase front counter
  else
    l          = size(Fr{nF},2);    % Number of solutions in front (nF)th
    idx        = front(nF).index;
    phie       = phi(:,idx);
    dis        = zeros(1,l);
    for j = 1:n_f
      [tmp1, idx1]  = sort(phie(j,:),'ascend');
      dis(idx1(1)) = Inf;
      dis(idx1(l)) = Inf;
      f_max        = tmp1(l);
      f_min        = tmp1(1);
      for k = 2:l-1
        temp1        = (tmp1(k+1) - tmp1(k-1))/(f_max - f_min);
        dis(idx1(k)) = dis(idx1(k)) + temp1;
      end
    end
  end
end

if ~isempty(dis)
  [tmp, idx2] = sort(dis,'descend');
  for i = kk:mu
    min_x(:,i)   = P(:,idx(idx2(i-kk+1)));
    min_sigma{i} = sigma_or{idx(idx2(i-kk+1))};
    min_alpha{i} = alpha_or{idx(idx2(i-kk+1))};
  end
end

end
%% END
