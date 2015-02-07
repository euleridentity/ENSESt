function dom_mat  = domination_matrix(n_viol, sum_viol, phi, m, n_f)
%%  dom_mat  = domination_matrix(n_viol, sum_viol, phi, m, n_f)
%   This function computes the 'Domination matrix' which specifies the
%   domination relation between two individuals using constrained-domination.
%
%   INPUT DATA:
%
%   - n_viol:   Vector which specifies the number of constraints violated
%               by each individual (1 x (mu+lambda) vector)
%   - sum_viol: Vector with the sum of the absolute values of the
%               constraints violated by each individual
%               (1 x (mu+lambda) vector)
%   - phi:      Vector matrix with the functions evaluated at each
%               population's individual.
%   - m:        Number of individuals (Positive integer number)
%   - n_f:      Number of functions to minimize (Positive integer number)
%
%   OUTPUT DATA:
%
%   - dom_mat: Domination matrix ((mu+lambda) x (mu+lambda) matrix)
%       dom_mat(i,j) =  1  : i dominates j
%       dom_mat(i,j) = -1  : j dominates i
%       dom_mat(p,q) =  0  : non-domination
%
%  BIBLIOGRAPHY:
%
%  - DEB, Kalyanmoy. "Multi-Objective optimization using evolutionary
%    algorithms". John Wiley & Sons, LTD. Kanpur, India. 2004.


%% Beginning
dom_mat  = zeros(m, m);

for i = 1:m-1
  for j = i+1:m
    %% 1. 'i' and 'j' are both feasible
    if (n_viol(i) == 0) && (n_viol(j)==0)
      i_dom_j = false;
      j_dom_i = false;
        for k = 1:n_f
          if phi(k,i) < phi(k,j)
            i_dom_j = true;
          elseif phi(k,i) > phi(k,j)
            j_dom_i = true;
          end
        end
        
        if i_dom_j && ~j_dom_i
          dom_mat(i,j) = 1;
        elseif ~i_dom_j && j_dom_i
          dom_mat(i,j) = -1;
        end

    %% 2. 'i' is feasible, and 'j' is infeasible
    elseif(n_viol(i) == 0) && (n_viol(j)~=0)
      dom_mat(i,j) = 1;

    %% 3. 'j' is feasible, and 'i' is infeasible
    elseif(n_viol(i)~= 0) && (n_viol(j)==0)
      dom_mat(i,j) = -1;

    %% 4. 'i' and 'j' are both infeasible
    else
      if (sum_viol(i) < sum_viol(j))
        dom_mat(i,j) = 1;
      elseif (sum_viol(i) > sum_viol(j))
        dom_mat(i,j) = -1;
      end
    end
  end
end

dom_mat = dom_mat - dom_mat';

end
%% END
