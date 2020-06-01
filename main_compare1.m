clear all
N=30;          %ï¿½ï¿½Ó¢ï¿½ï¿½ï¿½ï¿½Èºï¿½ï¿½Ä£
c1=1.5962;     %ï¿½ï¿½Ó¢ï¿½ï¿½ï¿½ï¿½Ñ§Ï°ï¿½ï¿½ï¿½ï¿½c1
c2=1.4962;     %ï¿½ï¿½Ó¢ï¿½ï¿½ï¿½ï¿½Ñ§Ï°ï¿½ï¿½ï¿½ï¿½c2
w1=0.89;       %ï¿½ï¿½Ó¢ï¿½ï¿½ï¿½Ó¹ï¿½ï¿½ï¿½È¨ï¿½ï¿½w
P=30;          %ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½Èºï¿½ï¿½Ä£
c11=1.4962;    %ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½Ñ§Ï°ï¿½ï¿½ï¿½ï¿½c1
c12=1.5962;    %ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½Ñ§Ï°ï¿½ï¿½ï¿½ï¿½c2
w11=0.7111;    %ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½Ó¹ï¿½ï¿½ï¿½È¨ï¿½ï¿½w      
D=30;          %ï¿½ï¿½ï¿½ï¿½ï¿½Õ¼ï¿½Î¬ï¿½ï¿½(ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½Êºï¿½3Î¬ï¿½ï¿½ï¿½ï¿½ï¿½Ï£ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿?1,2Î¬),ï¿½ï¿½ï¿½ï¿½Èºï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½10000*D
func_num=1;
runs=10;
I_strategy=1;
MaxDT=40000;    %Á£×ÓÈº×î´óµü´ú´ÎÊý
Max_FES=80000;
pg=zeros(MaxDT,6*28);
Pg=zeros(MaxDT,6,28);
fbest=zeros(runs,6*28);
Fbest=zeros(runs,6,28);
k=zeros(28,4);
f_mean=zeros(28);
f_std=zeros(28);
fhd1=str2func('benchmark_func1');
fhd=str2func('benchmark_func');
fhd2=str2func('benchmark_func_org');
for Ni=1:28
    func_num=Ni; 
     [gbest_trend_0,gbestval_0] = DE_DFPSO(fhd,N,P,c11,c12,w1,D,func_num,I_strategy,MaxDT,runs,Max_FES);
     [gbest_trend_1,gbestval_1,t1]= deopt(fhd,N+P,0.85,0.9,D,MaxDT,I_strategy,func_num,runs,Max_FES);
     [gbest_trend_2,gbestval_2,t2]= EPSO(fhd2,N+P,D,MaxDT,Max_FES,func_num,runs);
     [gbest_trend_3,gbestval_3,t3]= FDR(fhd,N+P,D,MaxDT,Max_FES,func_num,runs);
     [gbest_trend_4,gbestval_4]= CoDE(fhd1,N+P,D,runs,Max_FES,func_num);
%      [gbest_trend_4,gbestval_4,t4]= HPSO_TVAC(fhd,N+P,D,MaxDT,Max_FES,func_num);
     [gbest_trend_5,gbestval_5]=CODE_DFPSO(fhd,N,P,c11,c12,w1,D,func_num,runs,MaxDT,Max_FES);
%     [gbest_trend_5,gbestval_5,t5]= DMS_PSO_HS_func(fhd,200000,N,D,-100,100,MaxDT,func_num); %(feval,Max_FES,group_num,Particle_Number,Dimension,VRmin,VRmax,func_num)     
      for j=1:6
           s1='gbest_trend_';
           s2=num2str(j-1);
           s=strcat(s1,s2);
           Pg(:,j,Ni)=eval(s)';
           s1='gbestval_';
           s=strcat(s1,s2);
           Fbest(:,j,Ni)=eval(s)'; 
      end
end
for k=1:28
    str1='Sheet';
    str2=num2str(k);
    str=strcat(str1,str2);
    xlswrite('fbest_trend.xlsx',Pg(:,:,k),str) 
    xlswrite('fbest.xlsx',Pg(:,:,k),str) 
end
% xx=linspace(1,MaxDT,MaxDT);
% for k=1:8
%     figure
%     yy=Pg(:,xx,k);
%     plot(xx,yy,'--')
% end