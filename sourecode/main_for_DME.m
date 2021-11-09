clear;
clc;
set_num=5; % the number of dataset in the disease model
Re=zeros(6,18);
%% dataset
MAF=['0.05';'0.10';'0.20';'0.40'];
folder = '.\';
for p=1:6
    switch(p)
        case 1
            parameter='model = Additive; sample=1000; dimension = 4, H2 = 0.02, MAF = 0.05/0.10/0.20/0.40';
            path='.\DME_data\Additive\sample=1000\4_';
            dataFile = strcat(folder,'additive_4_1000.xls');
            H2 = '0.02';
        case 2
            parameter='model = Additive; sample=4000; dimension = 4, H2 = 0.02, MAF = 0.05/0.10/0.20/0.40';
            path='.\DME_data\Additive\sample=4000\4_';
            dataFile = strcat(folder,'additive_4_4000.xls');
            H2 = '0.02';
        case 3
            parameter='model = multiplicative; sample=1000; dimension = 4, H2 = 0.02, MAF = 0.05/0.10/0.20/0.40';
            path='.\DME_data\Multiplicative\sample=1000\4_';
            dataFile = strcat(folder,'multiplicative_4_1000.xls');
            H2 = '0.10';
        case 4
            parameter='model = multiplicative; sample=4000; dimension = 4, H2 = 0.02, MAF = 0.05/0.10/0.20/0.40';
            path='.\DME_data\Multiplicative\sample=4000\4_';
            dataFile = strcat(folder,'multiplicative_4_4000.xls');
            H2 = '0.10';
        case 5
            parameter='model = threshold; sample=1000; dimension = 4, H2 = 0.02, MAF = 0.05/0.10/0.20/0.40';
            path='.\DME_data\Threshold\sample=1000\4_';
            dataFile = strcat(folder,'threshold_4_1000.xls');
            H2 = '0.10';
        case 6
            parameter='model = threshold; sample=4000; dimension = 4, H2 = 0.02, MAF = 0.05/0.10/0.20/0.40';
            path='.\DME_data\Threshold\sample=4000\4_';
            dataFile = strcat(folder,'threshold_4_4000.xls');
            H2 = '0.10';
    end
        PARAMS.snp=100;
        PARAMS.dim_epi=4;
        PARAMS.dim_task=PARAMS.dim_epi+1;  %dimenson of task,while you need to set dim_epi<=dim_task
        PARAMS.num_ant=50;
        PARAMS.Max_FEs=2000*factorial(PARAMS.dim_epi)*PARAMS.dim_task;
        PARAMS.Esize=10;

    p_value=0.05/nchoosek(PARAMS.snp,PARAMS.dim_epi)*PARAMS.dim_task;
    Header={'Model','Power1','Power2', 'FEs','Total_FEs','RunTime','Total_time','TPR','SPC', 'PPV', 'ACC','FDR','F1'};
    sheet = 1;
    xlRange = 'A1';
    xlswrite(dataFile,Header,sheet,xlRange)    
    f=1;
    for m=1:1
        TP = 0;
        FP = 0;
        TN = 0;
        FN = 0;
        SFEs=0;
        ave_FEs=0;
        Runtime=0;
        Power1=0;
        Power2=0;
        Total_time=0;
        Model=strcat(MAF(m,:),'_',H2);
        filename=strcat(path,Model,'\');
        for i=1:set_num
            if i<10
                STR=['00',num2str(i),'.txt'];
            elseif i<100
                STR=['0',num2str(i),'.txt'];
            else
                STR=[num2str(i),'.txt'];
            end
            Filename=strcat(filename,STR);
            if strcmp(Filename(end-2:end),'txt')
                data=dlmread(Filename,'\t',1,0);
            elseif strcmp(Filename(end-2:end),'mat')
                load(Filename);
            end
            aim_snp=PARAMS.snp-PARAMS.dim_epi+1:1:PARAMS.snp; % functional SNP combination;
            [candidate,FEs,time,Total_FEs]=MTAco(data,PARAMS,aim_snp);
            if ismember(aim_snp,candidate,'rows')
                Power1=Power1+1;
            end
            G_set=[];
            for g=1:size(candidate,1)  
                G_score=Gtest_score(data(:,candidate(g,:)),data(:,end));
                if G_score < p_value
                    G_set = [G_set;candidate(g,:)];
                end
            end       
            TN = TN + size(candidate,1);
            if isempty(G_set)
                FN = FN + 1;      
            elseif ismember(aim_snp,G_set,'rows') 
                G_setSize =  length(G_set(:,1));
                Power2 = Power2 + 1;
                TP = TP + 1;
                TN = TN - G_setSize;
                fprintf('dataSet:%d    [    *****Success*****    Power1/Power2:%d/%d     FES=%3d     Total_FEs=%3d     time=%.3f]\n',i,Power1,Power2,FEs,Total_FEs,time);
            else
                G_setSize =  length(G_set(:,1));
                FP = FP + G_setSize;
                TN = TN - G_setSize;
                fprintf('dataSet:%d    [     *****Fail*****     Power1/Power2:%d/%d     FES=%3d     Total_FEs=%3d     time=%.3f]\n',i,Power1,Power2,FEs,Total_FEs,time);
            end
            
            SFEs=SFEs+FEs;
            Runtime=Runtime+time*(FEs/Total_FEs);
            Total_time = Total_time + time;
            ave_FEs = ave_FEs+Total_FEs;
        end
        SFEs=SFEs/set_num;
        ave_FEs=ave_FEs/set_num;
        Runtime=Runtime/set_num;
        Total_time=Total_time/set_num;
        TPR = TP/(TP + FN + 1e-10);
        SPC = TN/(FP + TN + 1e-10);
        PPV = TP/(TP + FP+ 1e-10);
        ACC = (TP + TN)/(TP + TN + FN + FP);
        FDR = 1 - PPV;
        F1 = 2*TP/(2*TP + FP + FN);
        Model=cellstr(Model);
        f=f+1;
        xlRange = strcat('A', num2str(f));
        xlswrite(dataFile,Model,sheet,xlRange)
        Results=[Power1,Power2,SFEs,ave_FEs, Runtime,Total_time,TPR,SPC,PPV,ACC,FDR,F1];
        xlRange = strcat('B', num2str(f));
        xlswrite(dataFile,Results,sheet,xlRange)
    end
end





