function [newPop, newPorp] = growCancer(pop, r, k, coeff)
%growCancer(population, r, k, coeff)   Calculates and returns new 
%population of cancer cells, using the G function.

%Vectorized verion. Error on the scale of 10^-7?
porp = pop./sum(pop);
E = coeff * porp';      %Calculates values of E, is 3x1
%mu = m/(deadliness+ r*v);

G = r .* (1 - sum(pop)*(1-E')./k);


%Unvectorized Version
% E = [0 0 0];
% G = [0 0 0];
% porp = pop./sum(pop);
% for i = 1:1:3
% 
%     E(i) = coeff(i,:) * porp';
% 
%     G(i) = r(i) * (k(i) - (1-E(i)) * sum(pop))/k(i);
% 
% end


%Stop G from being insane?
if (any(G == -inf))
    G(G==-inf) = -1;
end

%keep all pops alive???
if(any(pop<0.1))
    pop(pop<0.1) = 0.1;
    %disp('Pop Saved')
end

newPop = pop + pop .* G; %Find new populations
newPorp = newPop./sum(newPop);

%Debugging

% if pop ~= newPop
%     disp('AHHHH!!!')
%     disp([pop newPop])
%     disp(abs(pop-newPop))
%     pause(0.05)
%     %E(-1) = -1;
% else
%     disp('ok')
% end
end


