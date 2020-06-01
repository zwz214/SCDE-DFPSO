function [position,value,iteration]= EPSO(benchmark_func_org,num_particle,dimension,max_iteration,max_FES,func_num,runs,range)


% Description 
% CLPSO/Sa-PSO(PSO,FDR,HPSO,LIPS,CLPSO-gbest)

% CLPSO: c1=3~1.5, w1=0.9~0.4
% Sa-PSO: PSO/FDR/HPSO/LIPS/CLPSO-gbest, w=0.9~0.2
% PSO: c2_1=2.5~0.5, c2_2=0.5~2.5
% HPSO: c3_1=2.5~0.5, c3_2=0.5~2.5
% LIPS: nsize=3
% CLPSO with gbest:c4_1=2.5~0.5, c4_2=0.5~2.5 & Pc=0.5
% clear all;
% clc;
% max_FES=100000;
% max_iteration=5000;
% num_particle=20;
% dimension=10;
% func_num=3;
% range=[-100,100];
num_g=num_particle;
num_g1=8;
num_g2=num_particle-num_g1;
Pg=zeros(1,max_FES); 
y=(-1400:100:1400);
% load fbias_data
% y=f_bias;
outcome=[];
% result_data=zeros(1,max_FES);
% result=zeros(,num_g); % Initial fitness values
for time=1:runs
    rand('state',sum(100*clock));
    range_min=range(1)*ones(num_g,dimension); % Range for initial swarm's elements
    range_max=range(2)*ones(num_g,dimension);
    interval = range_max-range_min;
    v_max=interval*0.5; % 20% of interval
    v_min=-v_max;
    pos = range_min+ interval.*rand(num_g,dimension);    % Initial position
    vel =v_min+(v_max-v_min).*rand(num_g,dimension); % Initial velocity
    fitcount=0;
    fitold=1;
    % for i=1:num_g
    result=benchmark_func_org(pos',func_num);
    fitcount=fitcount+num_g;
    %      result_data=result(i);
    %      result_data(fitcount)=result(i);
    % end

    [gbest_val,g_index]=min(result);
    gbest_pos=pos(g_index,:); % Initial global best position

    pbest_pos=pos; % Initial pbest position
    pbest_val=result';

    % Method: CLPSO
    c1=3-(1:max_iteration)*1.5/max_iteration;
    w1=0.9-(1:max_iteration)*(0.5/max_iteration);

    % Sa-PSO (PSO/FDR/HPSO/LIPS/CLPSO-gbest)

    w2=0.9-(1:max_iteration)*(0.7/max_iteration);

    % Method 1: PSO
    c2_1=2.5-(1:max_iteration)*2/max_iteration;
    c2_2=0.5+(1:max_iteration)*2/max_iteration;

    % Method 2: FDR_PSO
    fii=[1 1 2];

    % Method 3: HPSO_TVAC
    c3_1=2.5-(1:max_iteration)*2/max_iteration;
    c3_2=0.5+(1:max_iteration)*2/max_iteration;
    re_init_vel=range(2)-(1:max_iteration)*(0.9*range(2))/max_iteration;


    % Method 4: LIPS
    nsize=3;

    % Method 5: CLPSO
    c4_1=2.5-(1:max_iteration)*2/max_iteration;
    c4_2=0.5+(1:max_iteration)*2/max_iteration;
    obj_func_slope=zeros(num_g,1);
    fri_best=(1:num_g)'*ones(1,dimension);
    j=0:(1/(num_g-1)):1; 
    j=j*10;
    Pc=ones(dimension,1)*(0.0+((0.5).*(exp(j)-exp(j(1)))./(exp(j(num_g))-exp(j(1))))); 

    for i=1:num_g1
        fri_best(i,:)=i*ones(1,dimension);
        friend1=ceil(num_g1*rand(1,dimension));
        friend2=ceil(num_g1*rand(1,dimension));
        friend=(pbest_val(friend1)<pbest_val(friend2)).*friend1'+(pbest_val(friend1)>=pbest_val(friend2)).*friend2';
        toss=ceil(rand(1,dimension)-Pc(:,i)');
        if toss==ones(1,dimension)
           temp_index=randperm(dimension);
           toss(1,temp_index(1))=0;
           clear temp_index;
        end
        fri_best(i,:)=(1-toss).*friend'+toss.*fri_best(i,:);
        for d=1:dimension
            fri_best_pos(i,d)=pbest_pos(fri_best(i,d),d);
        end
    end

    for i=num_g1+1:num_g1+num_g2
        fri_best(i,:)=i*ones(1,dimension);
        friend1=ceil(num_g*rand(1,dimension));
        friend2=ceil(num_g*rand(1,dimension));
        friend=(pbest_val(friend1)<pbest_val(friend2)).*friend1'+(pbest_val(friend1)>=pbest_val(friend2)).*friend2';
        toss=ceil(rand(1,dimension)-Pc(:,i)');
        if toss==ones(1,dimension)
           temp_index=randperm(dimension);
           toss(1,temp_index(1))=0;
           clear temp_index;
        end
        fri_best(i,:)=(1-toss).*friend'+toss.*fri_best(i,:);
        for d=1:dimension
            fri_best_pos(i,d)=pbest_pos(fri_best(i,d),d);
        end
    end


    k=0; % iteration 
    num_strategy = 5;
    LP=50; % Learning Period
    % temp_pos=[];

    while fitcount<=max_FES

            k=k+1; 
            gbest_pos_temp=repmat(gbest_pos,num_g,1);

    if fitcount<=(0.9*max_FES)

          % Group 1: CLPSO
          for i=1:num_g1 
                if obj_func_slope(i)>5
                           fri_best(i,:)=i*ones(1,dimension);
                           friend1=(ceil(num_g1*rand(1,dimension)));
                           friend2=(ceil(num_g1*rand(1,dimension)));
                           friend=(pbest_val(friend1)<pbest_val(friend2)).*friend1'+(pbest_val(friend1)>=pbest_val(friend2)).*friend2';
                           toss=ceil(rand(1,dimension)-Pc(:,i)');
                           if toss==ones(1,dimension)
                              temp_index=randperm(dimension);
                              toss(1,temp_index(1))=0;
                              clear temp_index;
                           end
                            fri_best(i,:)=(1-toss).*friend'+toss.*fri_best(i,:);
                            for d=1:dimension
                                fri_best_pos(i,d)=pbest_pos(fri_best(i,d),d);
                            end
                            obj_func_slope(i)=0;
                 end % if obj_func_slope(i)>5    

                delta(i,:)=(c1(1).*rand(1,dimension).*(fri_best_pos(i,:)-pos(i,:)));
                vel(i,:)=w1(k)*vel(i,:)+delta(i,:);
                vel(i,:)=((vel(i,:)<v_min(i,:)).*v_min(i,:))+((vel(i,:)>v_max(i,:)).*v_max(i,:))+(((vel(i,:)<v_max(i,:))&(vel(i,:)>v_min(i,:))).*vel(i,:));
                pos(i,:)=pos(i,:)+vel(i,:);

                if (sum(pos(i,:)>range_max(i,:))+sum(pos(i,:)<range_min(i,:))==0)
                    result(i)=benchmark_func_org(pos(i,:)',func_num);% Evaluate fitness   
                    fitcount=fitcount+1;
    %                  result_data=result(i);
    %                 result_data(fitcount)=result(i);
                    if fitcount>=max_FES
                        break;
                    end
                    if  result(i)<pbest_val(i) % update pbest value and position
                        pbest_pos(i,:)=pos(i,:);   
                        pbest_val(i)=result(i);
                        obj_func_slope(i)=0;
                    else
                        obj_func_slope(i)=obj_func_slope(i)+1;
                    end
                    if  pbest_val(i)<gbest_val; % update gbest value and postion
                        gbest_pos=pbest_pos(i,:); 
                        gbest_val=pbest_val(i);
                    end   
                end      
         end % CLPSO        

          for i=num_g1+1:num_g1+num_g2

                if k <=1
                    pk = ones(1,5)*1/num_strategy;
                    rk = 0:1/num_strategy:1;
                    success_mem = zeros(1,num_strategy);
                    failure_mem = zeros(1,num_strategy);
                    sk = zeros(1,num_strategy);
                elseif mod(fitcount,LP)==0
                    total = (success_mem+failure_mem);
                    total(find(total==0))=1;
                    sk = (success_mem./total)+0.01;
                    %sk_index = isnan(sk);
                    %sk(find(sk_index==1)) = 0.0001;
                    pk = sk./sum(sk);
    %                 fprintf('success = %f %f %f %f %f\n',success_mem);
    %                 fprintf('fail = %f %f %f %f %f\n',failure_mem);
    %                 fprintf('probability = %f %f %f %f %f %d\n',pk,k);
                    rk = [0 cumsum(pk)];
                    success_mem = zeros(1,num_strategy);
                    failure_mem = zeros(1,num_strategy);   
                end  
                probability = rand(1);
    %             prob_value(k)=probability;

                if probability>= rk(1) && probability < rk(2)
                    % PSO  
    %                 fprintf('PSO\n');
                    strategy_k = 1;
                    delta(i,:)=c2_1(k).*rand(1,dimension).*(pbest_pos(i,:)-pos(i,:))+c2_2(k).*rand(1,dimension).*(gbest_pos_temp(i,:)-pos(i,:));
                    vel(i,:)=w2(k).*vel(i,:)+delta(i,:); 

                elseif probability>= rk(2) && probability < rk(3)   
                    % FDR-PSO
    %                 fprintf('FDR_PSO\n')
                    strategy_k=2;
                    dis=abs(repmat(pbest_pos(i,:),num_g2,1)-pbest_pos(1:num_g2,:));
                    fiterr=repmat(pbest_val(i),1,num_g2)-pbest_val(1:num_g2)';
                    fiterr=repmat(fiterr',1,dimension);
                    fiterr=fiterr-(dis==zeros(num_g2,dimension)).*fiterr;
                    dis=dis+(dis==zeros(num_g2,dimension));
                    FDR=fiterr./dis;
                    [fdr,Fid]=max(FDR);
                    for dimcnt=1:dimension
                        Pnd(i,dimcnt)=pbest_pos(Fid(dimcnt),dimcnt);
                    end
                    delta(i,:)=fii(1).*rand(1,dimension).*(pbest_pos(i,:)-pos(i,:))+fii(2).*rand(1,dimension).*(gbest_pos_temp(i,:)-pos(i,:))+fii(3).*rand(1,dimension).*(Pnd(i,:)-pos(i,:));
                    vel(i,:)=w2(k).*vel(i,:)+delta(i,:); 


                elseif probability>= rk(3) && probability < rk(4)   
                    % HPSO
    %                 fprintf('HPSO\n');
                    strategy_k = 3;
                    vel(i,:)=(c3_1(k).*rand(1,dimension).*(pbest_pos(i,:)-pos(i,:)))+(c3_2(k).*rand(1,dimension).*(gbest_pos_temp(i,:)-pos(i,:)));
                    for d=1:dimension
                        if vel(i,d)==0
                            if (rand(1)<0.5)
                                vel(i,d)=rand(1)*re_init_vel(k);
                            else
                                vel(i,d)=-rand(1)*re_init_vel(k);
                            end
                        end
                        vel(i,d)=sign(vel(i,d))*min(abs(vel(i,d)),range(2));
                    end

                elseif probability>= rk(4) && probability < rk(5)
                        %LIPS
    %                     fprintf('LIPS\n');
                        strategy_k = 4;
                        EU_dist=dist(pos(i,:),pbest_pos'); 
                        EU_dist(i)=max(EU_dist);
                        [min_dist,min_index]=sort(EU_dist); 
                        fi=(4.1./nsize).*rand(nsize,dimension);
                        FIP=sum(fi.*pbest_pos(min_index(1:nsize),:))./sum(fi);  
                        delta(i,:)=sum(fi).*(FIP-pos(i,:));
                        vel(i,:)=0.7298.*(vel(i,:)+delta(i,:));
    %              
               elseif probability>= rk(5) && probability < rk(6)
                        %CLPSO
    %                     fprintf('CLPSO\n');
                        strategy_k = 5;                          
                        delta(i,:)=(c4_1(k).*rand(1,dimension).*(fri_best_pos(i,:)-pos(i,:)))+(c4_2(k).*rand(1,dimension).*(gbest_pos_temp(i,:)-pos(i,:)));
                        vel(i,:)=w2(k)*vel(i,:)+delta(i,:);
    %                     vel(i,:)=((vel(i,:)<v_min(i,:)).*v_min(i,:))+((vel(i,:)>v_max(i,:)).*v_max(i,:))+(((vel(i,:)<v_max(i,:))&(vel(i,:)>v_min(i,:))).*vel(i,:));
                         if obj_func_slope(i)>5
                           fri_best(i,:)=i*ones(1,dimension);
                           friend1=(ceil(num_g*rand(1,dimension)));
                           friend2=(ceil(num_g*rand(1,dimension)));
                           friend=(pbest_val(friend1)<pbest_val(friend2)).*friend1'+(pbest_val(friend1)>=pbest_val(friend2)).*friend2';
                           toss=ceil(rand(1,dimension)-Pc(:,i)');
                           if toss==ones(1,dimension)
                              temp_index=randperm(dimension);
                              toss(1,temp_index(1))=0;
                              clear temp_index;
                           end
                            fri_best(i,:)=(1-toss).*friend'+toss.*fri_best(i,:);
                            for d=1:dimension
                                fri_best_pos(i,d)=pbest_pos(fri_best(i,d),d);
                            end
                            obj_func_slope(i)=0;
                        end                                                            
                end
    %             strategy(k)=strategy_k;
                % for all
                vel(i,:)=((vel(i,:)<v_min(i,:)).*v_min(i,:))+((vel(i,:)>v_max(i,:)).*v_max(i,:))+(((vel(i,:)<v_max(i,:))&(vel(i,:)>v_min(i,:))).*vel(i,:));
                pos(i,:)=pos(i,:)+vel(i,:); 

                if (sum(pos(i,:)>range_max(i,:))+sum(pos(i,:)<range_min(i,:))==0)
                    result(i)=benchmark_func_org(pos(i,:)',func_num);% Evaluate fitness   
                    fitcount=fitcount+1;
    %                 result_data(fitcount)=result(i);
                    if fitcount>=max_FES
                        break;
                    end
                    if  result(i)<pbest_val(i) % update pbest value and position
                        pbest_pos(i,:)=pos(i,:);   
                        pbest_val(i)=result(i);
                        success_mem(strategy_k) = success_mem(strategy_k) +1;
                    else
                        failure_mem(strategy_k) = failure_mem(strategy_k) + 1;
                    end
                    if strategy_k==5 && result(i)<pbest_val(i)
                        obj_func_slope(i)=0;
                    else
                        obj_func_slope(i)=obj_func_slope(i)+1;
                    end
                    if  pbest_val(i)<gbest_val; % update gbest value and postion
                        gbest_pos=pbest_pos(i,:); 
                        gbest_val=pbest_val(i);
                    end   
                end      
          end

    else

         for i=1:num_g

                if k <=1
                    pk = ones(1,5)*1/num_strategy;
                    rk = 0:1/num_strategy:1;
                    success_mem = zeros(1,num_strategy);
                    failure_mem = zeros(1,num_strategy);
                    sk = zeros(1,num_strategy);
                elseif mod(fitcount,LP)==0
                    total = (success_mem+failure_mem);
                    total(find(total==0))=1;
                    sk = (success_mem./total)+0.01;
                    %sk_index = isnan(sk);
                    %sk(find(sk_index==1)) = 0.0001;
                    pk = sk./sum(sk);
    %                 fprintf('success = %f %f %f %f %f\n',success_mem);
    %                 fprintf('fail = %f %f %f %f %f\n',failure_mem);
    %                 fprintf('probability = %f %f %f %f %f %d\n',pk,k);
                    rk = [0 cumsum(pk)];
                    success_mem = zeros(1,num_strategy);
                    failure_mem = zeros(1,num_strategy);   
                end   
    %             pk_1(k)=pk(1);
    %             pk_2(k)=pk(2);
    %             pk_3(k)=pk(3);
    %             pk_4(k)=pk(4);
    %             pk_5(k)=pk(5);
                probability = rand(1);
    %             prob_value(k)=probability;

                if probability>= rk(1) && probability < rk(2)
                    % PSO  
    %                 fprintf('FDR_PSO\n');
                    strategy_k = 1;
                    delta(i,:)=c2_1(k).*rand(1,dimension).*(pbest_pos(i,:)-pos(i,:))+c2_2(k).*rand(1,dimension).*(gbest_pos_temp(i,:)-pos(i,:));
                    vel(i,:)=w2(k).*vel(i,:)+delta(i,:);              

                elseif probability>= rk(2) && probability < rk(3)
                    % FDR_PSO  
    %                 fprintf('FDR_PSO\n');
                    strategy_k = 2;
                    dis=abs(repmat(pbest_pos(i,:),num_g2,1)-pbest_pos(1:num_g2,:));
                    fiterr=repmat(pbest_val(i),1,num_g2)-pbest_val(1:num_g2)';
                    fiterr=repmat(fiterr',1,dimension);
                    fiterr=fiterr-(dis==zeros(num_g2,dimension)).*fiterr;
                    dis=dis+(dis==zeros(num_g2,dimension));
                    FDR=fiterr./dis;
                    [fdr,Fid]=max(FDR);
                    for dimcnt=1:dimension
                        Pnd(i,dimcnt)=pbest_pos(Fid(dimcnt),dimcnt);
                    end
                    delta(i,:)=fii(1).*rand(1,dimension).*(pbest_pos(i,:)-pos(i,:))+fii(2).*rand(1,dimension).*(gbest_pos_temp(i,:)-pos(i,:))+fii(3).*rand(1,dimension).*(Pnd(i,:)-pos(i,:));
                    vel(i,:)=w2(k).*vel(i,:)+delta(i,:); 

                elseif probability>= rk(3) && probability < rk(4)   
                    % HPSO
    %                 fprintf('HPSO\n');
                    strategy_k = 3;
                    vel(i,:)=(c3_1(k).*rand(1,dimension).*(pbest_pos(i,:)-pos(i,:)))+(c3_2(k).*rand(1,dimension).*(gbest_pos_temp(i,:)-pos(i,:)));
                    for d=1:dimension
                        if vel(i,d)==0
                            if (rand(1)<0.5)
                                vel(i,d)=rand(1)*re_init_vel(k);
                            else
                                vel(i,d)=-rand(1)*re_init_vel(k);
                            end
                        end
                        vel(i,d)=sign(vel(i,d))*min(abs(vel(i,d)),range(2));
                    end

                elseif probability>= rk(4) && probability < rk(5)
                        %LIPS
    %                     fprintf('FIPS\n');
                        strategy_k = 4;
                        EU_dist=dist(pos(i,:),pbest_pos'); 
                        EU_dist(i)=max(EU_dist);
                        [min_dist,min_index]=sort(EU_dist); 
                        fi=(4.1./nsize).*rand(nsize,dimension);
                        FIP=sum(fi.*pbest_pos(min_index(1:nsize),:))./sum(fi);  
                        delta(i,:)=sum(fi).*(FIP-pos(i,:));
                        vel(i,:)=0.7298.*(vel(i,:)+delta(i,:));
    %              
               elseif probability>= rk(5) && probability < rk(6)
                        %CLPSO
    %                     fprintf('CLPSO\n');
                        strategy_k = 5; 
                        delta(i,:)=(c4_1(k).*rand(1,dimension).*(fri_best_pos(i,:)-pos(i,:)))+(c4_2(k).*rand(1,dimension).*(gbest_pos_temp(i,:)-pos(i,:)));
                        vel(i,:)=w2(k)*vel(i,:)+delta(i,:);
    %                     vel(i,:)=((vel(i,:)<v_min(i,:)).*v_min(i,:))+((vel(i,:)>v_max(i,:)).*v_max(i,:))+(((vel(i,:)<v_max(i,:))&(vel(i,:)>v_min(i,:))).*vel(i,:));

                        if obj_func_slope(i)>5
                           fri_best(i,:)=i*ones(1,dimension);
                           friend1=(ceil(num_g*rand(1,dimension)));
                           friend2=(ceil(num_g*rand(1,dimension)));
                           friend=(pbest_val(friend1)<pbest_val(friend2)).*friend1'+(pbest_val(friend1)>=pbest_val(friend2)).*friend2';
                           toss=ceil(rand(1,dimension)-Pc(:,i)');
                           if toss==ones(1,dimension)
                              temp_index=randperm(dimension);
                              toss(1,temp_index(1))=0;
                              clear temp_index;
                           end
                            fri_best(i,:)=(1-toss).*friend'+toss.*fri_best(i,:);
                            for d=1:dimension
                                fri_best_pos(i,d)=pbest_pos(fri_best(i,d),d);
                            end
                            obj_func_slope(i)=0;
                        end                                      
                end
    %             strategy(k)=strategy_k;
                % for all
                vel(i,:)=((vel(i,:)<v_min(i,:)).*v_min(i,:))+((vel(i,:)>v_max(i,:)).*v_max(i,:))+(((vel(i,:)<v_max(i,:))&(vel(i,:)>v_min(i,:))).*vel(i,:));
                pos(i,:)=pos(i,:)+vel(i,:); 

                if (sum(pos(i,:)>range_max(i,:))+sum(pos(i,:)<range_min(i,:))==0)
                    result(i)=benchmark_func_org(pos(i,:)',func_num); % Evaluate fitness   
                    fitcount=fitcount+1;
    %                  result_data=result(i);
    %                 result_data(fitcount)=result(i);
                    if fitcount>=max_FES
                        break;
                    end
                    if  result(i)<pbest_val(i) % update pbest value and position
                        pbest_pos(i,:)=pos(i,:);   
                        pbest_val(i)=result(i);
                        success_mem(strategy_k) = success_mem(strategy_k) +1;
                    else
                        failure_mem(strategy_k) = failure_mem(strategy_k) + 1;
                    end
                    if strategy_k==5 && result(i)<pbest_val(i)
                        obj_func_slope(i)=0;
                    else
                        obj_func_slope(i)=obj_func_slope(i)+1;
                    end
                    if  pbest_val(i)<gbest_val % update gbest value and postion
                        gbest_pos=pbest_pos(i,:); 
                        gbest_val=pbest_val(i);
                    end   
                end
         end
    end
    if time==1&&fitcount<=max_FES
        Pg(1,fitold:fitcount)=gbest_val-y(func_num); 
    end
    fitold=fitcount;
       if fitcount>=max_FES
          break;
       end 

       if (k==max_iteration)&&(fitcount<max_FES)
           k=k-1;
       end

    end
    outcome=[outcome gbest_val-y(func_num)];
end
position=Pg;
value=outcome;
iteration=k;
% end