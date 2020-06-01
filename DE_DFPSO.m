%*************************************************************************%
%***********************粒子群动态跟随差分进化算法*******************************%
%*************************粒子向量定义规则*********************************%
%    1:D    D+1:2D     2D+1:3D     3D+1           3D+2           3D+3        D+4    %
%    位置    速度      历史最优    适应度   变异标志位/从众阈值   微变化次数    标志位  %
function [xm,fv] = DE_DFPSO(fitness,N,P,c11,c12,w11,D,func_num,I_strategy,MaxDT,runs,Max_FES,range)
format long;  
% flg=1;
%eps=10e-8;    %满足条件限制的误差(在不知道最小值时候不用设置)
MAXindependent=50;
b=0.6;    
c13=1.396;
F_weight=0.85;
o=zeros(1,MaxDT); 
Pg=zeros(1,Max_FES);    %迭代过程中的最优适应度值
F_CR=0.9;
y=(-1400:100:1400);
% load fbias_data
% y=f_bias;
% range=[-100,100];
outcome=[];
for TIME=1:runs
   rand('state',sum(100*clock));
    %初始化粒子的速度和位置
    % for i=1:N
    %     for j=1:2*D
    %       elite(i,j)=rand;%初始化速度
    %     end
    % end
    % elite=[(elite(:,1:D)-0.5)*200,elite(:,D+1:3*D+4)];
    fitcount=0;
    fitold=1;
    FVr_rot  = (0:1:N-1);      % rotating index array (size N)
    elite=zeros(N,3*D+4); %精英粒子信息矩阵
    public=zeros(P,3*D+4);%大众粒子信息矩阵，3*D+2列从众阈值
    range_min=range(1)*ones(N,D);
    range_max=range(2)*ones(N,D);
    interval=range_max-range_min;
    v_max=interval*0.5;
    v_min=-v_max;
    elite(:,1:D)=range_min+ interval.*rand(N,D);
    elite(:,D+1:2*D)=v_min+(v_max-v_min).*rand(N,D);
    elite(:,2*D+1:3*D)=elite(:,1:D);
    public(:,1:D)=range_min+ interval.*rand(N,D);
    public(:,D+1:2*D)=v_min+(v_max-v_min).*rand(N,D);
    public(:,2*D+1:3*D)=public(:,1:D);
    %初始化局部最优位置与粒子初始位置相同
    % for i=1:N
    %     for j=2*D+1:3*D
    %         elite(i,j)=elite(i,j-2*D);
    %     end
    % end
    % for i=1:P
    %     for j=1:2*D
    % 	  public(i,j)=rand;
    %     end
    % end
    % public=[(public(:,1:D)-0.5)*200,public(:,D+1:3*D+4)];
    %public=[public(:,1:D),public(:,D+1:3*D+4)];
    %初始化局部最优位置与粒子初始位置相同
    % for i=1:P
    %     for j=2*D+1:3*D
    % 		public(i,j)=public(i,j-2*D);
    %     end
    % end
    %计算精英粒子的适应度
    elite(:,3*D+1)=fitness(elite(:,2*D+1:3*D)',func_num);
    fitcount=fitcount+N;
    elite(:,3*D+2)=elite(:,3*D+1);
    %计算大众粒子的适应度，并为各大众粒子分配从众阈值
    public(:,3*D+1)=fitness(public(:,2*D+1:3*D)',func_num);
    fitcount=fitcount+P;
    public(:,3*D+3)=0;    %连续不变次数初始化为0
    public(:,3*D+4)=0;    %标志位初始化为0

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
    %初始化寄存器与全局最优位置
    NEXTelite=sortrows(elite,3*D+1);
    Ebest=NEXTelite(1,2*D+1:3*D);
    % Ebest_valold=Ebest_val;
    NEXTpublic=sortrows(public,3*D+1);
    Pbest=NEXTpublic(1,2*D+1:3*D);
    if NEXTelite(1,3*D+1)<NEXTpublic(1,3*D+1)
        Gbest=Ebest;
        Gbest_val=NEXTelite(1,3*D+1);
    else
        Gbest=Pbest;
        Gbest_val=NEXTpublic(1,3*D+1);
    end
    %%%%%%%%%进入主循环%%%%%%%%%%
    for time=1:MaxDT
        W=(w11-0.61)*(MaxDT-time)/MaxDT+0.61;
        pop=elite(:,1:D);
        popold=pop;
        ind=randperm(4);
        o(time)=N;
          FVr_a1  = randperm(N);                   % shuffle locations of vectors
          FVr_rt  = rem(FVr_rot+ind(1),N);         % rotate indices by ind(1) positions
          FVr_a2  = FVr_a1(FVr_rt+1);              % rotate vector locations
          FVr_rt  = rem(FVr_rot+ind(2),N);
          FVr_a3  = FVr_a2(FVr_rt+1);                
          FVr_rt  = rem(FVr_rot+ind(3),N);
          FVr_a4  = FVr_a3(FVr_rt+1);               
          FVr_rt  = rem(FVr_rot+ind(4),N);
          FVr_a5  = FVr_a4(FVr_rt+1);

          FM_bm = repmat(Gbest,o(time),1);       
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
                 if (rand < 0.5)                               % Pmu = 0.5
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
         fitcount=fitcount+1;
         if S<elite(k,3*D+2)
             elite(k,1:D)=FM_ui(k,:);
             elite(k,3*D+2)=S;
         end
         if  elite(k,3*D+2)<elite(k,3*D+1)
             elite(k,2*D+1:3*D)=elite(k,1:D);
             elite(k,3*D+1)=elite(k,3*D+2);
         end     
    end
    k=0;
    %更新全局最优粒子
    Ebest1=Ebest;
    elite=sortrows(elite,3*D+1);
    Ebest=elite(1,2*D+1:3*D);
    Ebest_val=elite(1,3*D+1);
       if Ebest_val<Gbest_val
            Gbest=Ebest;
            Gbest_val=Ebest_val;
       end

    %初始化影响半径
        gbest=repmat(Gbest,30,1);
        dist=sum(abs(gbest-public(:,1:D)));
        R=b*max(dist);
    %计算受全局最优粒子的影响的粒子
        m=0;
        for i=1:P
             if norm(NEXTpublic(i,1:D)-Gbest)<R  %&&NEXTpublic(i,3*D+3)>-MAXindependent)%||(norm(NEXTpublic(i,1:D)-Gbest)>R&&NEXTpublic(i,3*D+3)>MAXindependent)
    %          if norm(public(i,1:D)-Gbest)<=R
                public(i,3*D+4)=1;
                m=m+1;
            else
                public(i,3*D+4)=0; 
            end
        end
        public=sortrows(public,-3*D-4);
    %按从众阈值对范围内的粒子进行排序
        NEXTpublic=[sortrows(public(1:m,:),3*D+2);public(m+1:P,:)];
        sum1=NEXTpublic(1,3*D+2);
    %根据从众阈值挑选出跟随粒子，剩余的都归为独立粒子
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
    %从众粒子速度更新、位置更新、适应度更新和粒子最优位置更新
        for i=1:k
                public(i,D+1:2*D)=W*NEXTpublic(i,D+1:2*D)+c11*rand*(NEXTpublic(i,2*D+1:3*D)-NEXTpublic(i,1:D))+c12*rand*(Ebest(1:D)-NEXTpublic(i,1:D));
                public(i,1:D)=NEXTpublic(i,1:D)+public(i,1+D:2*D);
       %连续成为跟随粒子的次数
%            if public(i,3*D+3)<=0
%               public(i,3*D+3)=public(i,3*D+3)-1;
%            else
%                public(i,3*D+3)=0;
%            end
           %适应度更新
          if(sum(public(i,1:D)>range_max(i,:))+sum(public(i,1:D)<range_min(i,:))==0)
                newfitness=fitness(public(i,1:D)',func_num);
                     fitcount=fitcount+1;
                if newfitness<public(i,3*D+1)
                        public(i,2*D+1:3*D)=public(i,1:D);
                        public(i,3*D+1)=newfitness;
                end
          else
              m= public(i,1:D) ;
              m(m>=100)=range(1)+(range(2)-range(1))*rand;
              m(m<=-100)=range(1)+(range(2)-range(1))*rand;
              public(i,1:D)=m;
          end
        end
    %独立粒子更新
        for i=k+1:P
           for j=D+1:2*D
                 public(i,j)=W*NEXTpublic(i,j)+c11*rand*(NEXTpublic(i,j+D)-NEXTpublic(i,j-D))+c12*rand*(Ebest1(j-D)-NEXTpublic(i,j-D))+c13*rand*(Pbest(j-D)-NEXTpublic(i,j-D));
           end 
           for j=1:D
               public(i,j)=NEXTpublic(i,j)+public(i,j+D);
          %判断粒子是否飞出边界
            if public(i,j)>100
               public(i,j)=100;
            end
           end
         %连续成为独立粒子次数
%            if public(i,3*D+3)>=0
%               public(i,3*D+3)=public(i,3*D+3)+1;
%            else
%                public(i,3*D+3)=0;
%            end
            if(sum(public(i,1:D)>range_max(i,:))+sum(public(i,1:D)<range_min(i,:))==0)
               newfitness=fitness(public(i,1:D)',func_num);%适应度更新
                    fitcount=fitcount+1;
                if newfitness<public(i,3*D+1)
                        public(i,2*D+1:3*D)=public(i,1:D);
                        public(i,3*D+1)=newfitness;
                end
            else
              m= public(i,1:D) ;
              m(m>=100)=range(1)+(range(2)-range(1))*rand;
              m(m<=-100)=range(1)+(range(2)-range(1))*rand;
              public(i,1:D)=m;
            end
        end
     %更新全局最优粒子
    public=sortrows(public,3*D+1);
    Pbest=public(1,2*D+1:3*D);
       if public(1,3*D+1)<Gbest_val
            Gbest=Pbest;
            Gbest_val=public(1,3*D+1);
       end
       if TIME==1&&fitcount<=Max_FES
        Pg(1,fitold:fitcount)=Gbest_val-y(func_num);
       end
       fitold=fitcount;
        if fitcount>=Max_FES
            break;
        end
    %     if fitness(Gbest',func_num)-y(func_num)<eps
    %         flg=0;
    %         break;
    %     else
    %         flg=1;
    %     end
    end
    outcome=[outcome Gbest_val-y(func_num)];
end
%绘制进化代数和适应度关系曲线图
% xx=linspace(1,MaxDT,MaxDT);
xm=Pg;
fv=outcome;
% if flg==1
%     fv=Gbest_val-y(func_num);
% else
%     fv=0;
% end


    
   
