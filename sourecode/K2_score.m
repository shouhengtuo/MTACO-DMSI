function K2Score= K2_score(snp_com,state)
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

for i=1:xrow   %Count the frequency of each genotype combination
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

A = sum((Hcase - Hcontrol).^2);
sample = disease + control;
%% K2 score 
    K2score=0;

    for i=1:lrow
        if sample(i)>0
            y=My_factorial(sample(i)+1);
            r=My_factorial(disease(i))+My_factorial(control(i));
            K2score=K2score+(r-y);
        end
    end
   K2Score =abs(K2score);
   %% f is function used to calculate log form factorial
    function f=My_factorial(e)
        f=0;
        if e>0
            for o=1:e
                f=f+log(o);
            end
        end
    end
end
