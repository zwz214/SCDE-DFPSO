%*************************************************************************%
%***********************����Ӧ�������Ⱥ�㷨�����ò�ֽ������죩*******************************%
%*************************���������������*********************************%
%    1:D    D+1:2D     2D+1:3D     3D+1           3D+2           3D+3        D+4    %
%    λ��    �ٶ�      ��ʷ����    ��Ӧ��   �����־λ/������ֵ   ΢�仯����    ��־λ  %
function [xm,fv] = DFDEPSO(fitness,N,w1,P,c11,c12,w11,D,func_num,I_strategy,MaxDT)
format long;  
% %��ʼ����������
%k=1;           %�������Ӹ���������
% status=zeros(1,MaxDT);
% q=1;           %�������Ӹ�����������ӣ�1��Gnbest(1,:),2:Gnbest(2,:)
% E1best=Gnbest(1,:);   %
% elite1=zeros(N,3*D+4);%%%%%%%%%%%%%%%%����
% d=9;           %��ֵ�仯��������
% x=0.06;        %�����ʵ��ڲ���
% u=0.2;         %�����ʵ��ڲ���
flg=1;
%eps=10e-8;    %�����������Ƶ����(�ڲ�֪����Сֵʱ��������)
MAXindependent=50;
var=2;         %����Ȩֵ
b=0.6;         
Re=0;          %����ֵ���������´���
c3=0.55;
c1=2.5-(1:MaxDT)*2/MaxDT;
c2=0.5+(1:MaxDT)*2/MaxDT;
c13=1.396;
F_weight=0.85;
thr=zeros(1,MaxDT); 
Threshold=zeros(1,MaxDT);
o=zeros(1,MaxDT); 
Pg=zeros(1,MaxDT);    %���������е�������Ӧ��ֵ
Gnbest=rands(2,D);    %����Ӧ����õ�����λ�ñ���
var_elite=zeros(2,D); %�������ӼĴ���
elite=zeros(N,3*D+4); %��Ӣ������Ϣ����
public=zeros(P,3*D+4);%����������Ϣ����3*D+2�д�����ֵ
F_CR=0.9;
y=(-1400:100:1400);
FVr_rot  = (0:1:N-1);      % rotating index array (size N)
rand('state',sum(100*clock));
%��ʼ�����ӵ��ٶȺ�λ��
for i=1:N
    for j=1:2*D
      elite(i,j)=rand;%��ʼ���ٶ�
    end
end
elite=[(elite(:,1:D)-0.5)*200,elite(:,D+1:3*D+4)];
%elite=[elite(:,1:D),elite(:,D+1:3*D+4)];
%��ʼ���ֲ�����λ�������ӳ�ʼλ����ͬ
for i=1:N
    for j=2*D+1:3*D
        elite(i,j)=elite(i,j-2*D);
    end
end
for i=1:P
    for j=1:2*D
	  public(i,j)=rand;
    end
end
public=[(public(:,1:D)-0.5)*200,public(:,D+1:3*D+4)];
%public=[public(:,1:D),public(:,D+1:3*D+4)];
%��ʼ���ֲ�����λ�������ӳ�ʼλ����ͬ
for i=1:P
    for j=2*D+1:3*D
		public(i,j)=public(i,j-2*D);
    end
end
%���㾫Ӣ���ӵ���Ӧ��
elite(:,3*D+1)=fitness(elite(:,2*D+1:3*D)',func_num);
elite(:,3*D+2)=elite(:,3*D+1);
%����������ӵ���Ӧ�ȣ���Ϊ���������ӷ��������ֵ
public(:,3*D+1)=fitness(public(:,2*D+1:3*D)',func_num);
public(:,3*D+3)=0;    %�������������ʼ��Ϊ0
public(:,3*D+4)=0;    %��־λ��ʼ��Ϊ0
for i=1:P    
  if i<0.5*P
	    public(i,3*D+2)=1;
    elseif i<0.75*P
        public(i,3*D+2)=fix(P/3);
    elseif i<0.875*P
        public(i,3*D+2)=fix(0.5625*P);
    elseif i<0.9375*P  
        public(i,3*D+2)=fix(0.6563*P);
    else
        public(i,3*D+2)=fix(0.75*P);
  end
end
%��ʼ���Ĵ�����ȫ������λ��
NEXTelite=sortrows(elite,3*D+1);
Ebest=NEXTelite(1,2*D+1:3*D);
Ebest_val=NEXTelite(1,3*D+1);
Ebest_valold=Ebest_val;
NEXTpublic=sortrows(public,3*D+1);
Pbest=NEXTpublic(1,2*D+1:3*D);
if NEXTelite(3*D+1)<NEXTpublic(3*D+1)
    Gbest=Ebest;
else
    Gbest=Pbest;
end
%%%%%%%%%������ѭ��%%%%%%%%%%
for time=1:MaxDT
 k=1;
      if rem(time,100)<50
         q=1;
     else
         q=2;
     end  
%  Threshold(1,time)=(1-time/MaxDT)^d*(fitness(ones(1,D),D));
%  Threshold(1,time)=350*time^(-0.9185);
% Threshold(1,time)=1056*exp(-0.107*time);%D=10
% Threshold(1,time)=1706*exp(-0.06853*time);%D=30
%  Threshold(1,time)=2750*exp(-0.08565*time);%D=50
%Threshold(1,time)=7.842*exp(-0.03624*time);%fintness3
l=abs(fix((rand-0.5)*N))+1;
Threshold(time)=norm(NEXTelite(l,1:D)-Ebest);
 NEXTelite=elite;
 w=(w1-0.71)*(MaxDT-time)/MaxDT+0.71;
 W=(w11-0.61)*(MaxDT-time)/MaxDT+0.61;
 
%��Ӣ�����ٶȸ���
    for i=1:N
        for j=D+1:2*D
             elite(i,j)=w*NEXTelite(i,j)+c1(time)*rand*(NEXTelite(i,j+D)-NEXTelite(i,j-D))+c2(time)*rand*(Gbest(j-D)-NEXTelite(i,j-D))+c3*rand*(Ebest(j-D)-NEXTelite(i,j-D));
%             elite(i,j)=w*NEXTelite(i,j)+(c1+(Re/1000)^2)*rand*(NEXTelite(i,j+D)-NEXTelite(i,j-D))+(c2-(Re/1500)^2)*rand*(Gbest(j-D)-NEXTelite(i,j-D));
        end
    end
%��Ӣ����λ�ø��£�������������ʷ����λ��
    for i=1:N
        elite(i,1:D)=NEXTelite(i,1:D)+elite(i,1+D:2*D);
        for j=1:D
            if elite(i,j)>100 %�ж������Ƿ�ɳ��߽�
             elite(i,j)=100;
           end
        end
        elite(i,3*D+2)=fitness(elite(i,1:D)',func_num); %��Ӧ�ȸ���
        if   elite(i,3*D+2)<elite(i,3*D+1)
             elite(i,2*D+1:3*D)=elite(i,1:D); 
             elite(i,3*D+1)=elite(i,3*D+2);
        end
    end
% %����ȫ������λ��
E1best=Ebest;
elite=sortrows(elite,3*D+1);
Ebest=elite(1,2*D+1:3*D);
Ebest_val=elite(1,3*D+1);
   if Ebest_val<Ebest_valold
        Gbest=Ebest;
        Pg(time)=elite(1,3*D+1);
   end
NEXTelite=elite;
    %ͳ�����������´���
thr(1,time)=Ebest_valold-Ebest_val;
if abs(thr(1,time))==0&&Re<500
    Re=Re+1;
elseif abs(thr(1,time))==0&&Re>=500
    Re=500;
else
    Re=0;
end
%�ж��Ƿ������������  
% if (u+Re*x)>2*rand
if Re>5&&time>500
       for i=1:N
         if (norm(NEXTelite(i,1:D)-Ebest))<(Threshold(time))
           NEXTelite(i,3*D+2)=i;
           o(time)=o(time)+1;
%            o1=o(time);
        else
           NEXTelite(i,3*D+2)=0;            
        end
      end
%�����ӽ��н������
    NEXTelite=sortrows(NEXTelite,-3*D-2);
    elite=NEXTelite;
%******�����ֽ���������**********%
    pop=elite(:,1:D);
    popold=pop;
    ind=randperm(4);
      FVr_a1  = randperm(N);                   % shuffle locations of vectors
      FVr_rt  = rem(FVr_rot+ind(1),N);         % rotate indices by ind(1) positions
      FVr_a2  = FVr_a1(FVr_rt+1);              % rotate vector locations
      FVr_rt  = rem(FVr_rot+ind(2),N);
      FVr_a3  = FVr_a2(FVr_rt+1);                
      FVr_rt  = rem(FVr_rot+ind(3),N);
      FVr_a4  = FVr_a3(FVr_rt+1);               
      FVr_rt  = rem(FVr_rot+ind(4),N);
      FVr_a5  = FVr_a4(FVr_rt+1);
 
      FM_bm = repmat(Ebest,o(time),1);       
      FM_pm1 = popold(FVr_a1(1:o(time)),:);             % shuffled population 1
      FM_pm2 = popold(FVr_a2(1:o(time)),:);             % shuffled population 2
      FM_pm3 = popold(FVr_a3(1:o(time)),:);             % shuffled population 3
      FM_pm4 = popold(FVr_a4(1:o(time)),:);             % shuffled population 4
      FM_pm5 = popold(FVr_a5(1:o(time)),:);             % shuffled population 5
          
      FM_mui = rand(o(time),D) < F_CR;  % all random numbers < F_CR are 1, 0 otherwise
      FM_mpo = FM_mui < 0.5;
        if (I_strategy == 1)                             % DE/rand/1
            FM_ui = FM_pm3 + F_weight*(FM_pm1 - FM_pm2);   % differential variation
            FM_ui = popold(1:o(time),:).*FM_mpo + FM_ui.*FM_mui;     % crossover
            FM_origin = FM_pm3;
        elseif (I_strategy == 2)                         % DE/local-to-best/1
            FM_ui = popold(1:o(time),:) + F_weight*(FM_bm-popold(1:o(time),:)) + F_weight*(FM_pm1 - FM_pm2);
            FM_ui = popold(1:o(time),:).*FM_mpo + FM_ui.*FM_mui;
            FM_origin = popold(1:o(time),:);
        elseif (I_strategy == 3)                         % DE/best/1 with jitter
            FM_ui = FM_bm + (FM_pm1 - FM_pm2).*((1-0.9999)*rand(o(time),D)+F_weight);               
            FM_ui = popold(1:o(time),:).*FM_mpo + FM_ui.*FM_mui;
            FM_origin = FM_bm;
        elseif (I_strategy == 4)                         % DE/rand/1 with per-vector-dither
             f1 = ((1-F_weight)*rand(N,1)+F_weight);
             for k=1:D
                FM_pm5(:,k)=f1;
             end
             FM_ui = FM_pm3 + (FM_pm1 - FM_pm2).*FM_pm5;    % differential variation
             FM_origin = FM_pm3;
             FM_ui = popold(1:o(time),:).*FM_mpo + FM_ui.*FM_mui;     % crossover
         elseif (I_strategy == 5)                          % DE/rand/1 with per-vector-dither
             f1 = ((1-F_weight)*rand+F_weight);
             FM_ui = FM_pm3 + (FM_pm1 - FM_pm2)*f1;         % differential variation
             FM_origin = FM_pm3;
             FM_ui = popold(1:o(time),:).*FM_mpo + FM_ui.*FM_mui;     % crossover
          else                                              % either-or-algorithm
             if (rand < 0.5);                               % Pmu = 0.5
                FM_ui = FM_pm3 + F_weight*(FM_pm1 - FM_pm2);% differential variation
             else                                           % use F-K-Rule: K = 0.5(F+1)
                FM_ui = FM_pm3 + 0.5*(F_weight+1.0)*(FM_pm1 + FM_pm2 - 2*FM_pm3);
             end
             FM_origin = FM_pm3;
             FM_ui = popold(1:o(time),:).*FM_mpo + FM_ui.*FM_mui;     % crossover     
        end
%*****************SELETION**********************%
for k=1:o(time)
     for j=1:D %----boundary constraints via bounce back-------
        if (FM_ui(k,j) > 100)
           FM_ui(k,j) = 100 + rand*(FM_origin(k,j) - 100);
        end
        if (FM_ui(k,j) < -100)
           FM_ui(k,j) = -100 + rand*(FM_origin(k,j) + 100);
        end   
     end
     S=fitness(FM_ui(k,:)',func_num);
     if S<elite(k,3*D+2)
         elite(k,1:D)=FM_ui(k,:);
         elite(k,3*D+2)=S;
     end
     if elite(k,3*D+2)<elite(k,3*D+1)
         elite(k,2*D+1:3*D)=elite(k,1:D);
         elite(k,3*D+1)=elite(k,3*D+2);
     end     
end
end
k=0;
%����ȫ����������
elite=sortrows(elite,3*D+1);
Ebest=elite(1,2*D+1:3*D);
Ebest_valold=elite(1,3*D+1);
%    if elite(1,3*D+1)<Pg(time-1)
%         Gbest=Ebest;
%    end
Gnbest(1,:)=Gbest;
Gnbest(2,:)=elite(1,2*D+1:3*D);
    %��ʼ��Ӱ��뾶
R=b*norm(Ebest-Pbest)+b*norm(Ebest+Pbest);
%������ȫ���������ӵ�Ӱ�������
    m=0;
    for i=1:P
         if (norm(NEXTpublic(i,1:D)-Gbest)<R&&NEXTpublic(i,3*D+3)>-MAXindependent)%||(norm(NEXTpublic(i,1:D)-Gbest)>R&&NEXTpublic(i,3*D+3)>MAXindependent)
%          if norm(public(i,1:D)-Gbest)<=R
            public(i,3*D+4)=1;
            m=m+1;
        else
            public(i,3*D+4)=0; 
        end
    end
    public=sortrows(public,-3*D-4);
%��������ֵ�Է�Χ�ڵ����ӽ�������
    NEXTpublic=[sortrows(public(1:m,:),3*D+2);public(m+1:P,:)];
    sum1=NEXTpublic(1,3*D+2);
%���ݴ�����ֵ��ѡ���������ӣ�ʣ��Ķ���Ϊ��������
    for i=2:m
        if NEXTpublic(i-1,3*D+2)==NEXTpublic(i,3*D+2)
           sum1=sum1+NEXTpublic(i,3*D+2);
           k=k+1;
        elseif sum1>NEXTpublic(i,3*D+2)
            sum1=sum1+NEXTpublic(i,3*D+2);
            k=k+1;
        else
            break;
        end
    end
public=NEXTpublic;
%���������ٶȸ��¡�λ�ø��¡���Ӧ�ȸ��º���������λ�ø���
    for i=1:k
            public(i,D+1:2*D)=W*NEXTpublic(i,D+1:2*D)+c11*rand*(NEXTpublic(i,2*D+1:3*D)-NEXTpublic(i,1:D))+c12*rand*(Ebest(1:D)-NEXTpublic(i,1:D));
%             public(i,D+1:2*D)=W*NEXTpublic(i,D+1:2*D)+(c11+(Re/1000)^2)*rand*(NEXTpublic(i,2*D+1:3*D)-NEXTpublic(i,1:D))+(c12-(Re/1500)^2)*rand*(Gbest(1:D)-NEXTpublic(i,1:D));
            public(i,1:D)=NEXTpublic(i,1:D)+public(i,1+D:2*D);
   %������Ϊ�������ӵĴ���
       if public(i,3*D+3)<=0
          public(i,3*D+3)=public(i,3*D+3)-1;
       else
           public(i,3*D+3)=0;
       end
       %��Ӧ�ȸ���
        newfitness=fitness(public(i,1:D)',func_num);
        if newfitness<public(i,3*D+1)
                public(i,2*D+1:3*D)=public(i,1:D);
                public(i,3*D+1)=newfitness;
        end
    end
%�������Ӹ���
    for i=k+1:P
       for j=D+1:2*D
             public(i,j)=W*NEXTpublic(i,j)+c11*rand*(NEXTpublic(i,j+D)-NEXTpublic(i,j-D))+c12*rand*(E1best(j-D)-NEXTpublic(i,j-D))+c13*rand*(Pbest(j-D)-NEXTpublic(i,j-D));
%             public(i,j)=W*NEXTpublic(i,j)+(c11+(Re/1000)^2)*rand*(NEXTpublic(i,j+D)-NEXTpublic(i,j-D))+(c12-(Re/1500)^2)*rand*(Gnbest(q,j-D)-NEXTpublic(i,j-D))+(c13+(Re/500)^2)*rand*(Pbest(j-D)-NEXTpublic(i,j-D));
       end 
       for j=1:D
           public(i,j)=NEXTpublic(i,j)+public(i,j+D);
      %�ж������Ƿ�ɳ��߽�
        if public(i,j)>100
           public(i,j)=100;
        end
       end
     %������Ϊ�������Ӵ���
       if public(i,3*D+3)>=0
          public(i,3*D+3)=public(i,3*D+3)+1;
       else
           public(i,3*D+3)=0;
       end
       newfitness=fitness(public(i,1:D)',func_num);%��Ӧ�ȸ���
        if newfitness<public(i,3*D+1)
                public(i,2*D+1:3*D)=public(i,1:D);
                public(i,3*D+1)=newfitness;
        end
    end
 %����ȫ����������
public=sortrows(public,3*D+1);
Pbest=public(1,2*D+1:3*D);
   if public(1,3*D+1)<fitness(Gbest',func_num)
        Gbest=Pbest;
   end
    Pg(time)=fitness(Gbest',func_num);
%     if fitness(Gbest',func_num)-y(func_num)<eps
%         flg=0;
%         break;
%     else
%         flg=1;
%     end
end
%���ƽ�����������Ӧ�ȹ�ϵ����ͼ
xx=linspace(1,MaxDT,MaxDT);
xm=Pg(xx)-ones(1,MaxDT)*y(func_num);
if flg==1
    fv=fitness(Gbest',func_num)-y(func_num);
else
    fv=0;
end


    
   
