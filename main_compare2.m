clear all
N=30;          %��Ӣ����Ⱥ��ģ
c1=1.5962;     %��Ӣ����ѧϰ����c1
c2=1.4962;     %��Ӣ����ѧϰ����c2
w1=0.89;       %��Ӣ���ӹ���Ȩ��w
P=30;          %��������Ⱥ��ģ
c11=1.4962;    %��������ѧϰ����c1
c12=1.5962;    %��������ѧϰ����c2
w11=0.7111;    %�������ӹ���Ȩ��w      
D=30;          %�����ռ�ά��(�������ʺ�3ά�����ϣ��������?1,2ά),����Ⱥ��������10000*D
func_num=1;
runs=10;
I_strategy=1;
MaxDT=40000;    %����Ⱥ�������� 
Max_FES=800000;
fhd1=str2func('benchmark_func1');
fhd=str2func('benchmark_func');
fhd2=str2func('benchmark_func_org');
 parfor Ni=1:25 
    func_num=Ni; 
    switch Ni
        case {1,2,3,4,5,6,12,14}
            range=[-100,100];
        case 7
            range=[0,600];
        case 8
            range=[-32,32]; 
        case {9,10,15,16,17,18,19,20,21,22,23,24}
            range=[-5,5];
        case 11
            range=[-0.5,0.5];
        case 13
            range=[-3,1];
        case 25
            range=[-2,5];
    end
     [gbest_trend_0,gbestval_0] = DE_DFPSO(fhd,N,P,c11,c12,w1,D,func_num,I_strategy,MaxDT,runs,Max_FES,range);
     [gbest_trend_1,gbestval_1,t1]= deopt(fhd,N+P,0.85,0.9,D,MaxDT,I_strategy,func_num,runs,Max_FES,range);
     [gbest_trend_2,gbestval_2,t2]= EPSO(fhd2,N+P,D,MaxDT,Max_FES,func_num,runs,range);
     [gbest_trend_3,gbestval_3,t3]= FDR(fhd,N+P,D,MaxDT,Max_FES,func_num,runs,range);
     [gbest_trend_4,gbestval_4]= CoDE(fhd1,N+P,D,runs,MaxDT,Max_FES,func_num);
%      [gbest_trend_4,gbestval_4,t4]= HPSO_TVAC(fhd,N+P,D,MaxDT,Max_FES,func_num);
     [gbest_trend_5,gbestval_5]=CSDE_DFPSO(fhd,N,P,c11,c12,w1,D,func_num,runs,MaxDT,Max_FES,range);
%      [gbest_trend_6,gbestval_6]=CoDE_DFPSO(fhd,N,P,c11,c12,w1,D,func_num,runs,MaxDT,Max_FES,range);
%     [gbest_trend_5,gbestval_5,t5]= DMS_PSO_HS_func(fhd,200000,N,D,-100,100,MaxDT,func_num); %(feval,Max_FES,group_num,Particle_Number,Dimension,VRmin,VRmax,func_num)     
     data=[gbest_trend_0',gbest_trend_1',gbest_trend_2',gbest_trend_3',gbest_trend_4',gbest_trend_5'];
     data1=[gbestval_0',gbestval_1',gbestval_2',gbestval_3',gbestval_4',gbestval_5'];
%        data=gbest_trend_6';
%       data1=gbestval_6';
    str1='Sheet';
    str2=num2str(Ni);
    str=strcat(str1,str2);
    xlswrite('fbest_trend_CODE.xlsx',data,str) 
    xlswrite('fbest_val_CODE.xlsx',data1,str) 
end