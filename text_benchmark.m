clear all
N=20;          %��Ӣ����Ⱥ��ģ
c1=1.4962;     %��Ӣ����ѧϰ����c1
c2=1.4962;     %��Ӣ����ѧϰ����c2
w=0.7298;      %��Ӣ���ӹ���Ȩ��w
P=20;          %��������Ⱥ��ģ
c11=1.4962;    %��������ѧϰ����c1
c12=1.4962;    %��������ѧϰ����c2
w1=0.6111;     %�������ӹ���Ȩ��w      
D=30;          %�����ռ�ά��(�������ʺ�3ά�����ϣ��������1,2ά),����Ⱥ����������10000*D
MaxDT=5000;    %����Ⱥ����������
Max_FES=10000;
func_num=1;                    % I_strategy     1 --> DE/rand/1:
runs=10;                       %                2 --> DE/local-to-best/1
I_strategy=3;                  %                3 --> DE/best/1 with jitter
Pg=zeros(runs,MaxDT,28);      
fbest=zeros(28,runs);
f_mean=zeros(28);
f_std=zeros(28);
fhd=str2func('benchmark_func');
for Ni=1:6
    func_num=Ni;
    for Nj=1:runs
%             [gbest_trend,gbestval]=DE_DFPSO(fhd,N,P,c11,c12,w1,D,func_num,I_strategy,MaxDT,Max_FES);
           [gbest_trend_1,gbestval_1]=CODE_DFPSO(fhd,N,P,c11,c12,w1,D,func_num,MaxDT,Max_FES);
%            [gbest_trend_1,gbestval_1,t1]= deopt(fhd,N+P,0.85,0.9,D,MaxDT,I_strategy,func_num,Max_FES);
%           [gbest_trend_2,gbestval_2,t2]= EPSO(fhd,N+P,D,MaxDT,Max_FES,func_num);
%          [gbest_trend_3,gbestval_3,t3]= FDR(fhd,N+P,D,MaxDT,Max_FES,func_num);
%          [gbest_trend_4,gbestval_4,t4]= HPSO_TVAC(fhd,N+P,D,MaxDT,Max_FES,func_num);
        fbest(Ni,Nj)=gbestval;
        Pg(Nj,:,Ni)=gbest_trend;
        Pg(Nj+runs,:,Ni)=gbest_trend_1;
    end
%     f_mean(Ni)=mean(fbest(Ni,:));
%     f_std(Ni)=std(fbest(Ni,:));
end
xlswrite('fbest.xlsx',Pg(:,:,1)','Sheet1')
xlswrite('fbest.xlsx',Pg(:,:,2)','Sheet2')
xlswrite('fbest.xlsx',Pg(:,:,3)','Sheet3')
xlswrite('fbest.xlsx',Pg(:,:,4)','Sheet4')
xlswrite('fbest.xlsx',Pg(:,:,5)','Sheet5')
xlswrite('fbest.xlsx',fbest','Sheet6')
% xlswrite('fbest.xlsx',f_mean,'Sheet2')
% xlswrite('fbest.xlsx',f_std,'Sheet3')

