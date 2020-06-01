%*************************************************************************%
%***********************����Ӧ�������Ⱥ�㷨*******************************%
%*************************���������������*********************************%
%    1:D    D+1:2D     2D+1:3D     3D+1           3D+2           3D+3        D+4    %
%    λ��    �ٶ�      ��ʷ����    ��Ӧ��   �����־λ/������ֵ   ΢�仯����    ��־λ  %
function [xm,fv] = SASPSO(fitness,N,w1,P,c11,c12,w11,D,func_num,MaxDT)
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
thr=zeros(1,MaxDT); 
Threshold=zeros(1,MaxDT);
o=zeros(1,MaxDT); 
Pg=zeros(1,MaxDT);    %���������е�������Ӧ��ֵ
Gnbest=rands(2,D);    %����Ӧ����õ�����λ�ñ���
var_elite=zeros(2,D); %�������ӼĴ���
elite=zeros(N,3*D+4); %��Ӣ������Ϣ����
public=zeros(P,3*D+4);%����������Ϣ����3*D+2�д�����ֵ
y=(-1400:100:1400);
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
E1best=Ebest;
% %����ȫ������λ��
elite=sortrows(elite,3*D+1);
Ebest=elite(1,2*D+1:3*D);
   if elite(1,3*D+1)<fitness(Gbest',func_num)
        Gbest=Ebest;
   end
NEXTelite=elite;
    %ͳ�����������´���
thr(1,time)=fitness(E1best',func_num)-elite(1,3*D+1);
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
           o1=o(time);
        else
           NEXTelite(i,3*D+2)=0;            
        end
      end
%�����ӽ��н������
    NEXTelite=sortrows(NEXTelite,-3*D-2);
    elite=NEXTelite;
    for i=1:ceil(o(time)/2)
        e=rand;
        var_elite(1,1:D)=0.5*(NEXTelite(2*i-1,1:D)+NEXTelite(2*i-1,2*D+1:3*D))*e+0.5*(NEXTelite(2*i,1:D)+NEXTelite(2*i,2*D+1:3*D))*(1-e);
        var_elite(2,1:D)=0.5*(NEXTelite(2*i-1,1:D)+NEXTelite(2*i-1,2*D+1:3*D))*(1-e)+0.5*(NEXTelite(2*i,1:D)+NEXTelite(2*i,2*D+1:3*D))*e;
        fit1=fitness(var_elite(1,:)',func_num);
        fit2=fitness(var_elite(2,:)',func_num);
        if fit1<elite(2*i-1,3*D+1)
             elite(2*i-1,1:D)=var_elite(1,:);  
             elite(2*i-1,2*D+1:3*D)=var_elite(1,:);
             elite(2*i-1,3*D+1)=fit1;     
        end
        if fit2<elite(2*i,3*D+1)
            elite(2*i,1:D)=var_elite(2,:);
            elite(2*i,2*D+1:3*D)=var_elite(2,:);
            elite(2*i,3*D+1)=fit2;
        end
        if fit1>elite(2*i-1,3*D+1)&&fit2>elite(2*i,3*D+1)
            %����������������
               var_elite(1,:)=NEXTelite(2*i-1,1:D)+(1-time/MaxDT)^var*(fitness(100*ones(1,D)',func_num)-NEXTelite(2*i-1,1:D));
               var_elite(2,:)=NEXTelite(2*i-1,1:D)-(1-time/MaxDT)^var*(NEXTelite(2*i-1,1:D)-fitness(-100*ones(1,D)',func_num));
           if fitness(var_elite(1,:)',func_num)<fitness(var_elite(2,:)',func_num)
                  elite(2*i-1,1:D)=var_elite(1,:);
           else
                  elite(2*i-1,1:D)=var_elite(2,:);
           end      
               var_elite(1,:)=NEXTelite(2*i,1:D)+(1-time/MaxDT)^var*(fitness(100*ones(1,D)',func_num)-NEXTelite(2*i,1:D));
               var_elite(2,:)=NEXTelite(2*i,1:D)-(1-time/MaxDT)^var*(NEXTelite(2*i,1:D)-fitness(-100*ones(1,D)',func_num));            
            if fitness(var_elite(1,:)',func_num)<fitness(var_elite(2,:)',func_num)
                elite(2*i,1:D)=var_elite(1,:);
            else
                elite(2*i,1:D)=var_elite(2,:);
            end
        end
    end  
end
%����ȫ����������
elite=sortrows(elite,3*D+1);
Ebest=elite(1,2*D+1:3*D);
   if elite(1,3*D+1)<fitness(Gbest',func_num)
        Gbest=Ebest;
   end
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
            public(i,D+1:2*D)=W*NEXTpublic(i,D+1:2*D)+c11*rand*(NEXTpublic(i,2*D+1:3*D)-NEXTpublic(i,1:D))+c12*rand*(Gbest(1:D)-NEXTpublic(i,1:D));
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
             public(i,j)=W*NEXTpublic(i,j)+c11*rand*(NEXTpublic(i,j+D)-NEXTpublic(i,j-D))+c12*rand*(Gnbest(q,j-D)-NEXTpublic(i,j-D))+c13*rand*(Pbest(j-D)-NEXTpublic(i,j-D));
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


    
   
