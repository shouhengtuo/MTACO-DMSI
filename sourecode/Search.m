function pop=Search(PARAMS,pop)
% PARAMS
snp=PARAMS.snp;
dim_task=PARAMS.dim_task;
T0=0.8;
% Aco search
for dim=2:dim_task+1
    num_ant = size(pop(dim-1).SNPs,1);
    Randpos = randperm(snp);
    SNPs=ones(num_ant,dim);
    % Randomly initialize the starting position
    if num_ant<= snp
        for e = 1:num_ant
            SNPs(e,1) = Randpos(e);
        end
    else
        for e = 1:snp
            SNPs(e,1) = Randpos(e);
        end
        for e = snp+1:num_ant
            SNPs(e,1) = randperm(snp,1);
        end
    end

    Tau = pop(dim-1).Tau;
    for i=1:num_ant
        temp_Tau=Tau;
        temp_Tau(SNPs(i,1:dim))=0;
        for j = 2:dim
            p=temp_Tau;
            p=p/(sum(p));
            pcum = cumsum(p);
            mRate = rand;
            if mRate >= T0         % Roulette wheel selection
                select = ceil(rand*snp);
                while(ismember(select,SNPs(i,1:dim-1)))
                    select = ceil(rand*snp);
                end
            else                   
                select = find(pcum >= rand,1); % Random selection
            end
            SNPs(i,j) = select;
            temp_Tau(select)=0; % Set the pheromones of the selected SNP to 0
        end
    end
    [~,ia,~] = unique(sort(SNPs,2),'rows');
    pop(dim-1).SNPs=SNPs(ia,:);
end