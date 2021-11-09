function JS = JS_score(snp_com,state)
%%
 ua = unique(state);
if length(ua)~=2
    disp(' Class state is not equal to 2 !!!');
elseif min(ua)>0
    state = state -min(ua);
end

[xrow,xcol] = size(snp_com);
[Data,~,cid]=unique(snp_com,'rows');
[lrow,~]=size(Data);

Hcase = zeros(3,xcol);
Hcontrol = zeros(3,xcol);
sample=zeros(lrow,1);
disease=sample;
control=sample;

for i=1:xrow   %% Count the frequency of each genotype combination
    if state(i) == 1
        disease(cid(i)) = disease(cid(i)) + 1;
        for j = 1:xcol
            if snp_com(i,j) == 0
                Hcase(1,j) = Hcase(1,j) + 1;
            elseif snp_com(i,j) == 1
                Hcase(2,j) = Hcase(2,j) + 1;
            else
                Hcase(3,j) = Hcase(3,j) + 1;
            end
        end
    else
        control(cid(i)) = control(cid(i)) + 1;
        for j = 1:xcol
            if snp_com(i,j) == 0
                Hcontrol(1,j) = Hcontrol(1,j) + 1;
            elseif snp_com(i,j) == 1
                Hcontrol(2,j) = Hcontrol(2,j) + 1;
            else
                Hcontrol(3,j) = Hcontrol(3,j) + 1;
            end
        end
    end
    
end
%% Jensen-Shannon
Pdisease = disease / sum(disease);
Pcontrol = control / sum(control);
JS1 = 0;
JS2 = 0;
for i = 1:lrow
    if Pdisease(i)>0
        JS1 = JS1 + Pdisease(i).*log2(2*Pdisease(i) / (Pdisease(i) + Pcontrol(i)) );
    end
    if Pcontrol(i)>0
        JS2 = JS2 + Pcontrol(i).*log2(2*Pcontrol(i) / (Pdisease(i) + Pcontrol(i)) );
    end
end
% 0 <=JS1 + JS2<= lrow
 JS = 1/ (JS1 + JS2);