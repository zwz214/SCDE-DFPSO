
function [position,value,iteration]= HPSO_TVAC(benchmark,num_particle,dimension,max_iteration,max_FES,func_num)
format long;  

range=[-100,100];
Pg=zeros(1,max_iteration); 
y=(-1400:100:1400);
num_g=num_particle;

% result_data=zeros(1,max_FES);
result=zeros(1,num_particle); % Initial fitness values
% check_vel=[];


c1=2.5-(1:max_iteration)*2/max_iteration;
c2=0.5+(1:max_iteration)*2/max_iteration;

range_min=range(1)*ones(num_g,dimension); % Range for initial swarm's elements
range_max=range(2)*ones(num_g,dimension);
interval = range_max-range_min;
%v_max=interval*0.2; % 20% of interval
%v_min=-v_max;
re_init_vel=range(2)-(1:max_iteration)*(0.9*range(2))/max_iteration;

pos = range_min+interval.*rand(num_g,dimension);    % Initial position
%vel =v_min+(v_max-v_min).*rand(num_g,dimension); % Initial velocity


fitcount=0;

for i=1:num_g
     result(i)=benchmark(pos(i,:)',func_num);
     fitcount=fitcount+1;
%      result_data(fitcount)=result(i);
end


[gbest_val,g_index]=min(result);
gbest_pos=pos(g_index,:); % Initial global best position

pbest_pos=pos; % Initial pbest position
pbest_val=result;

k=0; % iteration 
% temp_pos=[];

while k<=max_iteration && fitcount<=max_FES
        
        k=k+1; 
        gbest_pos_temp=repmat(gbest_pos,num_g,1);
        vel=(c1(k).*rand(num_g,dimension).*(pbest_pos-pos))+(c2(k).*rand(num_g,dimension).*(gbest_pos_temp-pos));
        for i=1:num_g
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
         end
         pos=pos+vel;
        
        
        for i=1:num_g   % Evaluate Fitness
               result(i)=benchmark(pos(i,:)',func_num);% Evaluate fitness
               fitcount=fitcount+1;
%                result_data(fitcount)=result(i);
                if fitcount>=max_FES
                   break;
                end
            if  result(i)<pbest_val(i)   % update pbest value and position
                pbest_pos(i,:)=pos(i,:);   
                pbest_val(i)=result(i);
            end
            
            if  pbest_val(i)<gbest_val   % update gbest value and postion
                gbest_pos=pbest_pos(i,:); 
                gbest_val=pbest_val(i);
            end   
        end 
     Pg(k)=gbest_val-y(func_num);   
   if fitcount>=max_FES
      break;
   end 
        
   if (k==max_iteration)&&(fitcount<max_FES)
       k=k-1;
   end      
                
end
xx=linspace(1,max_iteration,max_iteration);
position=Pg(xx);
value=gbest_val-y(func_num);
iteration=k;
% max_FES=fitcount;
end


