function pop=updatePheromones2(pop,k)
%  update Pheromones
rou=0.8; % Pheromone Evaporation Factor 
Delta1=0.2;  % Pheromone increment from elite solutions
Delta2=0.1;  % Pheromone increment from Current population 
Theta = 0.1;
Tau = pop(k).Tau;
max_ESscore = max(pop(k).ES.score);
[min_ESscore,~] = min(pop(k).ES.score);
Tau= Tau.*(1-rou);
% Pheromone increment from Current population 
max_score = max(pop(k).score);
min_score = min(pop(k).score);
Ed = max_score-min_score+1e-10; %extremedifference
for i = 1:size(pop(k).SNPs,1)
    lambda1 = (max_score-pop(k).score(i,1))/ Ed; % normalized
    Tau(pop(k).SNPs(i,:),1) = Tau(pop(k).SNPs(i,:),1)+ lambda1*Delta2;
end
% Pheromone increment from elite solutions
for i = 1:size(pop(k).ES.Cand,1)
    Ed = max_ESscore-min_ESscore+1e-10; %extremedifference
    lambda2 = (max_ESscore-pop(k).ES.score(i,1))/ Ed; % normalized
    Tau(pop(k).ES.Cand(i,:),1) = Tau(pop(k).ES.Cand(i,:),1)+ lambda2*Delta1;
end
%
for i = 1:size(Tau)
    Tau(i,1) = max(Tau(i,1),Theta);
end
 pop(k).Tau = Tau;