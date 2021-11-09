clear;
clc;
set_num=5;
Re=zeros(6,18);
%%
dataFile = '.\NDME.xls';
for p=1:8
    switch(p)
        case 1
            path='.\NDME_data\NDME-1\Dim3_Id';
            PARAMS.dim_epi=3;
        case 2
            path='.\NDME_data\NDME-2\Dim3_Id';
            PARAMS.dim_epi=3;
        case 3
            path='.\NDME_dataNDME-3\Dim4_Id';
            PARAMS.dim_epi=4;
        case 4
            path='.\NDME_data\NDME-4\Dim4_Id';
            PARAMS.dim_epi=4;
        case 5
            path='.\NDME_data\NDME-5\Dim4_Id';
            PARAMS.dim_epi=4;
        case 6
            path='.\NDME_data\NDME-6\Dim5_Id';
            PARAMS.dim_epi=5;
        case 7
            path='.\NDME_data\NDME-7\Dim5_Id';
            PARAMS.dim_epi=5;
        case 8
            path='.\NDME_datas\NDME-8\Dim5_Id';
            PARAMS.dim_epi=5;
    end
    PARAMS.snp=100;
    PARAMS.num_ant=50;
    PARAMS.dim_task=2*PARAMS.dim_epi-1;  
    PARAMS.Max_FEs=60000*(PARAMS.dim_epi-1)^2;
    PARAMS.Esize=10;
    p_value=0.05/nchoosek(PARAMS.snp,PARAMS.dim_epi);
    Header={'Model','Power1','Power2', 'FEs','Total_FEs','RunTime','Total_time','TPR','SPC', 'PPV', 'ACC','FDR','F1'};
    sheet = 1;
    xlRange = 'A1';
    xlswrite(dataFile,Header,sheet,xlRange) 
    %% evaluate criterion
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
    for i=1:set_num
        Filename=strcat(path,num2str(i),'.mat');
        load(Filename);
        aim_snp=PARAMS.snp-PARAMS.dim_epi+1:1:PARAMS.snp;
        [candidate,FEs,time,Total_FEs]=MTAco(data,PARAMS,aim_snp);
        fprintf('dataSet:%d           FES=%3d     Total_FEs=%3d     time=%.3f]\n',i,FEs,Total_FEs,time);
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
    %%
    Model=strcat('NDME_',num2str(p));
    Model=cellstr(Model);
    xlRange = strcat('A', num2str(p+1));
    xlswrite(dataFile,Model,sheet,xlRange)
    Results=[Power1,Power2,SFEs,ave_FEs, Runtime,Total_time,TPR,SPC,PPV,ACC,FDR,F1];
    xlRange = strcat('B', num2str(p+1));
    xlswrite(dataFile,Results,sheet,xlRange);
    
end
    





