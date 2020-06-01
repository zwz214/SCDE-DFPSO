clear all
mex cec13_func.cpp -DWINDOWS
N=30;          %��Ӣ����Ⱥ��ģ
c1=1.5962;     %��Ӣ����ѧϰ����c1
c2=1.4962;     %��Ӣ����ѧϰ����c2
w1=0.89;       %��Ӣ���ӹ���Ȩ��w
P=30;          %��������Ⱥ��ģ
c11=1.4562;    %��������ѧϰ����c1
c12=1.5962;    %��������ѧϰ����c2
w11=0.7111;    %�������ӹ���Ȩ��w      
D=30;          %�����ռ�ά��(�������ʺ�3ά�����ϣ��������?1,2ά),����Ⱥ��������10000*D
func_num=1;
runs=10;
I_strategy=1;
MaxDT=40000;    %����Ⱥ����������
Max_FES=300000;
Pg =zeros(4,MaxDT,6);
fbest=zeros(28,4);
k=zeros(28,4);
range=[-100,100];
f_mean=zeros(28);
f_std=zeros(28);
fhd=str2func('cec13_func');
for Ni=1:28
        func_num=Ni; 
        [gbest_trend_0,gbestval_0] = DE_DFPSO(fhd,N,P,c11,c12,w1,D,func_num,I_strategy,MaxDT,runs,Max_FES,range);
        [gbest_trend_1,gbestval_1,t1]= deopt(fhd,N+P,0.85,0.9,D,MaxDT,I_strategy,func_num,runs,Max_FES,range);
        [gbest_trend_2,gbestval_2,t2]= EPSO(fhd,N+P,D,MaxDT,Max_FES,func_num,runs,range);
        [gbest_trend_3,gbestval_3,t3]= FDR(fhd,N+P,D,MaxDT,Max_FES,func_num,runs,range);
%       [gbest_trend_4,gbestval_4,t4]= HPSO_TVAC(fhd,N+P,D,MaxDT,Max_FES,func_num);
        [gbest_trend_4,gbestval_4]= CoDE(fhd,N+P,D,runs,MaxDT,Max_FES,func_num);
        [gbest_trend_5,gbestval_5]= SCDE_DFPSO1(fhd,N,P,c11,c12,w1,D,func_num,runs,MaxDT,Max_FES,range);
        data =[gbest_trend_0' gbest_trend_1' gbest_trend_2' gbest_trend_3' gbest_trend_4' gbest_trend_5'];
        data1 =[gbestval_0' gbestval_1' gbestval_2' gbestval_3' gbestval_4' gbestval_5'];
        str1='Sheet';
        str2=num2str(Ni);
        str=strcat(str1,str2);
        xlswrite('fbest_scde_cec2013_2.xlsx',data,str)
        xlswrite('fbest_scde_cec2013_val_2.xlsx',data1,str)
end

