function[P_value]=Gtest_score(snp_com,state)
ua = unique(state);
if length(ua)~=2 
    disp(' Class state is not equal to 2 !!!');
elseif min(ua)>0
    state = state -min(ua);
end
[xrow,xcol] = size(snp_com);
[Data,idx,cid]=unique(snp_com,'rows');
[lrow,~]=size(Data);
%% 

F=zeros(3,lrow+1);  %Column 1 holds the number of cases samples, Column 1 holds the number of controls samples, and column 3 holds the total number of samples.
E=zeros(2,lrow);  %% Expectation
% Count the frequency of each genotype combination
for i=1:xrow   
       if state(i)==0 %% case        
            F(1,cid(i))=F(1,cid(i))+1;        
       else 
            F(2,cid(i))=F(2,cid(i))+1;  
       end  
    
end



F(3,1:lrow)=sum(F(1:2,1:lrow),1);
F(:,lrow+1)=sum(F(:,1:lrow),2);


G=0;
Degree=(2-1)*(lrow-1);

for i=1:2
    for j=1:lrow
        O=F(i,j)/xrow;
        E=( ( F(i,lrow+1) * F(3,j) )/xrow )/xrow;
       
      if F(3,j)>2
          if O>0
              G=G+(xrow*O)*log(O/E);
          end
      
      elseif Degree>1
         Degree=Degree-0.5; 
      end
    end    
end
G=2*G;
P_value=1-chi2cdf(G,Degree);