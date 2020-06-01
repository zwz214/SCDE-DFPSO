function [gbesttrend,gbestval,iteration]= DMS_PSO_funcec08(feval,err0,Max_FES,Max_interation,group_num,Particle_Number,Dimension,Xmin,Xmax,VRmin,VRmax,func_num)
%usage:     [gbest,gbestval,fitcount,fit_cut,get_flag]= DMS_PSO_func08(1,0,0/-inf,500000,30,3,100,-100,100,-100,100)
%this program used to solve minimization problem, if you want wot solve a
%maximization problem please add "-" in your function.
%Outputs:
%           gbest:      the final result
%           gbestval:   the corresponding 
%           fitcount:   the corresponding used FEs
%           fit_cut:    the used FEs when achieve gbestval<err0
%           get_flag:   get_flag=1 means get the predefined criterion gbestval<err0
%Inputs:
%           fhd:        your problem function 
%           err0:       criterion, the FEs used for achieve this criterion will be recorded
%           err1:       stopping criterion  (please set err0/err1 to -Inf if you do not want the algorithm stop before getting Max_FEs)
%           Max_FEs:    The number of max function evaluations
%           group_num:  The group number of DMS-PSO
%           Particle_Number:    The number of particles in EACH group
%           Dimension:  the dimension of your problem
%           Xmin,Xmax:  the bounds of initilization
%           VRmin,VRmax:the bounds of search space
%           varargin:   all the parameters you need in your problem function

%javaclasspath('FractalFunctions.jar')
gbest_record=[];
rng('shuffle')
group_ps=Particle_Number;
ps=group_num*group_ps;
D=Dimension;
L_FES=3*D;
L_num=ceil(0.25.*group_num);
cc=[1.49445 1.49445];   %acceleration constants
iwt=0.729;
Pg=zeros(1,Max_interation); 
y=(-1400:100:1400);
% iwt=0.9-(1:me)*(0.7/me);

if length(VRmin)==1
    VRmin=repmat(VRmin,1,D);
    VRmax=repmat(VRmax,1,D);
end

VRmin=repmat(VRmin,ps,1);
VRmax=repmat(VRmax,ps,1);

if length(Xmin)==1
    Xmin=repmat(Xmin,1,D);
    Xmax=repmat(Xmax,1,D);
end
Xmin=repmat(Xmin,ps,1);
Xmax=repmat(Xmax,ps,1);

mv=0.2*(Xmax-Xmin);
Vmin=-mv;
Vmax=-Vmin;

pos=Xmin+(Xmax-Xmin).*rand(ps,D);
fhd='benchmark_func';
for i=1:ps;
    e(i,1)=feval(pos(i,:)',func_num);
end


fitcount=ps;
vel=Vmin+2.*Vmax.*rand(ps,D);%initialize the velocity of the particles
pbest=pos;
pbestval=e; %initialize the pbest and the pbest's fitness value
for i=1:group_num
    group_id(i,:)=[((i-1)*group_ps+1):i*group_ps];
    pos_group(group_id(i,:))=i;
    [gbestval(i),gbestid]=min(pbestval(group_id(i,:)));
    gbest(i,:)=pbest(group_id(i,gbestid),:);%initialize the gbest and the gbest's fitness value
end
% get_flag=0;
i=0;
while fitcount<0.95*Max_FES&&i<Max_interation
    
    i=i+1;
    for k=1:ps
    aa(k,:)=cc(1).*rand(1,D).*(pbest(k,:)-pos(k,:))+cc(2).*rand(1,D).*(gbest(pos_group(k),:)-pos(k,:));
    vel(k,:)=iwt.*vel(k,:)+aa(k,:); 
    vel(k,:)=(vel(k,:)>mv(k,:)).*mv(k,:)+(vel(k,:)<=mv(k,:)).*vel(k,:); 
    vel(k,:)=(vel(k,:)<(-mv(k,:))).*(-mv(k,:))+(vel(k,:)>=(-mv(k,:))).*vel(k,:);
    pos1(k,:)=pos(k,:)+vel(k,:); 
    keep_d=rand(1,D)<0.5;
    pos(k,:)=keep_d.*pbest(k,:)+(1-keep_d).*pos1(k,:); 
    if (sum(pos(k,:)>VRmax(k,:))+sum(pos(k,:)<VRmin(k,:)))==0;

    e(k,1)=feval(pos(k,:)',func_num);
    fitcount=fitcount+1;
    tmp=(pbestval(k)<e(k));
    temp=repmat(tmp,1,D);
    pbest(k,:)=temp.*pbest(k,:)+(1-temp).*pos(k,:);
    pbestval(k)=tmp.*pbestval(k)+(1-tmp).*e(k);%update the pbest
    if pbestval(k)<gbestval(pos_group(k))
    gbest(pos_group(k),:)=pbest(k,:);
    gbestval(pos_group(k))=pbestval(k);
    end
    end
    
    end
    
    change_flag=0;
%     for k=1:group_num
%         flag=1;
%         for kk=1:group_ps
%             flag=flag*(sum(abs(pbest(group_id(k,kk),:)-gbest(k,:))<1e-3)==D);
%         end
%         if flag==1
%             change_flag=1;break;
%         end
%     end
    
    
    if mod(i,100)==0
        options = optimset('LargeScale','off','MaxFunEvals',L_FES,'Display','off');
        [tmp,tmpid]=sort(gbestval);
        for k=1:L_num
            [x,fval,exitflag,output] = fminunc(fhd,gbest(tmpid(k),:),options,func_num);
            fitcount=fitcount+output.funcCount;
            if fval<gbestval(tmpid(k))
                [gbestval(tmpid(k)),gbestid]=min(pbestval(group_id(tmpid(k),:)));
                pbest(group_id(tmpid(k),gbestid),:)=x;
                pbestval(group_id(tmpid(k),gbestid))=fval;
                gbest(tmpid(k),:)=x;
                gbestval(tmpid(k))=fval;
            end
        end
    end
    
    if change_flag==1||mod(i,5)==0
        rc=randperm(ps);
        group_id=[];gbest=[];gbestval=[];
        for k=1:group_num
            group_id(k,:)=rc(((k-1)*group_ps+1):k*group_ps);
            pos_group(group_id(k,:))=k;
            [gbestval(k),gbestid]=min(pbestval(group_id(k,:)));
            gbest(k,:)=pbest(group_id(k,gbestid),:);
        end
    end
%     if round(i/100)==i/100
% %     plot(pos(:,D-1),pos(:,D),'b*');
% %     hold on
% %     plot(gbest(:,D-1),gbest(:,D),'r*');   
% %     hold off
% %     axis([VRmin(1,D-1),VRmax(1,D-1),VRmin(1,D),VRmax(1,D)])
% %     title(['PSO: ',num2str(i),' generations, Gbestval=',num2str(min(gbestval))]);  
% %     drawnow
% i,gbestval
%     end
% if min(gbestval)<err0&&get_flag==0
%     fit_cut=fitcount;
%     get_flag=1;
% end
% if min(gbestval)<err1
%     break;
% end

if mod(i,1)==0
    [aaa,bbb]=sort(gbestval);
    ggg=gbestval(bbb(1));
    gbest_record=[gbest_record;ggg];
end
aaa=[];bbb=[];ggg=[];
end
[tmp,tmpid]=sort(gbestval);
gbest=gbest(tmpid(1),:);
gbestval=gbestval(tmpid(1));

% if gbestval>err1
options = optimset('LargeScale','off','MaxFunEvals',0.05*Max_FES,'Display','off');
% options = optimset('LargeScale','off','MaxFunEvals',max_gen,'DiffMinChange',1e-20);
[x,fval,exitflag,output] = fminunc(fhd,gbest,options,func_num);
fitcount=fitcount+output.funcCount;
if fval<gbestval
    gbest=x;
    gbestval=fval;
end
% end
[aaa,bbb]=sort(gbestval);
ggg=gbestval(bbb(1));
gbest_record=[gbest_record;ggg];
aaa=[];bbb=[];ggg=[];

while fitcount<Max_FES
    
    i=i+1;
    for k=1:ps
    aa(k,:)=cc(1).*rand(1,D).*(pbest(k,:)-pos(k,:))+cc(2).*rand(1,D).*(gbest-pos(k,:));
    vel(k,:)=iwt.*vel(k,:)+aa(k,:); 
    vel(k,:)=(vel(k,:)>mv(k,:)).*mv(k,:)+(vel(k,:)<=mv(k,:)).*vel(k,:); 
    vel(k,:)=(vel(k,:)<(-mv(k,:))).*(-mv(k,:))+(vel(k,:)>=(-mv(k,:))).*vel(k,:);
    pos(k,:)=pos(k,:)+vel(k,:); 
    if (sum(pos(k,:)>VRmax(k,:))+sum(pos(k,:)<VRmin(k,:)))==0;
    e(k,1)=feval(pos(k,:)',func_num);
    fitcount=fitcount+1;
    tmp=(pbestval(k)<e(k));
    temp=repmat(tmp,1,D);
    pbest(k,:)=temp.*pbest(k,:)+(1-temp).*pos(k,:);
    pbestval(k)=tmp.*pbestval(k)+(1-tmp).*e(k);%update the pbest
    if pbestval(k)<gbestval
    gbest=pbest(k,:);
    gbestval=pbestval(k);
    end
    end
    
    end
    

% if gbestval<err0&&get_flag==0
%     fit_cut=fitcount;
%     get_flag=1;
% end

if mod(i,1)==0
    [aaa,bbb]=sort(gbestval);
    ggg=gbestval(bbb(1));
    gbest_record=[gbest_record;ggg];
end
aaa=[];bbb=[];ggg=[];

end


if gbestval<err0
    fit_cut=fitcount;
end
% 
% if get_flag==0
%     fit_cut=fitcount;
% end
