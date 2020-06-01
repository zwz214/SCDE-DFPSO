function [gbest,gbestval,fitcount]= DMS_PSO_HS_func(feval,Max_FES,group_num,Particle_Number,Dimension,VRmin,VRmax,func_num)

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

% gbest_record=[];
rand('state',sum(100*clock));
group_ps=Particle_Number;
ps=group_num*group_ps;
D=Dimension;
L_FES=100;
L_num=ceil(0.25.*group_num);
cc=[1.49445 1.49445];   %acceleration constants
iwt=0.729;
% iwt=0.9-(1:me)*(0.7/me);

if length(VRmin)==1
    VRmin=repmat(VRmin,1,D);
    VRmax=repmat(VRmax,1,D);
end
VRmin=repmat(VRmin,ps,1);
VRmax=repmat(VRmax,ps,1);
mv=0.2*(VRmax-VRmin);
Vmin=-mv;
Vmax=-Vmin;
pos=VRmin+(VRmax-VRmin).*rand(ps,D);

for i=1:ps;
    e(i,1)=feval(fhd,pos(i,:)',func_num);
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
get_flag=0;
i=0;
while fitcount<0.95*Max_FES
    
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
   
    if change_flag==1|mod(i,5)==0
        rc=randperm(ps);
        group_id=[];gbest=[];gbestval=[];
        for k=1:group_num
            group_id(k,:)=rc(((k-1)*group_ps+1):k*group_ps);
            pos_group(group_id(k,:))=k;
            [gbestval(k),gbestid]=min(pbestval(group_id(k,:)));
            gbest(k,:)=pbest(group_id(k,gbestid),:);
        end
    end

% Harmony search
if mod(i,1)==0
HMCR=0.98;
PAR=0.01+(0.99-0.01)/Max_FES*fitcount;
% PAR=0.2;
% PSOR=PAR;
for k=1:group_num
    for d=1:D
        if rand(1,1)<=HMCR
            rpd=randperm(group_ps);
            posHS(k,d)=pbest(group_id(k,rpd(1)),d);
            if rand(1,1)<=PAR
                posHS(k,d)=gbest(k,d);
            end
        else posHS(k,d)=VRmin(1,d)+(VRmax(1,d)-VRmin(1,d)).*rand(1,1);
        end
    end
    if (sum(posHS(k,:)>VRmax(k,:))+sum(posHS(k,:)<VRmin(k,:)))==0;
        eHS(k,1)=feval(posHS(k,:)',func_num);
        fitcount=fitcount+1;
        [gworstval(k),gworstid]=max(pbestval(group_id(k,:)));
        tmp=(gworstval(k)<eHS(k));
        temp=repmat(tmp,1,D);
        pbest(group_id(k,gworstid),:)=temp.*pbest(group_id(k,gworstid),:)+(1-temp).*posHS(k,:);
        pbestval(group_id(k,gworstid))=tmp.*pbestval(group_id(k,gworstid))+(1-tmp).*eHS(k);%update the pbest by HS
        if pbestval(group_id(k,gworstid),:)<gbestval(pos_group(k))
        gbest(pos_group(k),:)=pbest(group_id(k,gworstid),:);
        gbestval(pos_group(k))=pbestval(group_id(k,gworstid));
        end
    end
end
end

end
[tmp,tmpid]=sort(gbestval);
gbest=gbest(tmpid(1),:);
gbestval=gbestval(tmpid(1));

while fitcount<Max_FES
    
    i=i+1;
    for k=1:ps
    aa(k,:)=cc(1).*rand(1,D).*(pbest(k,:)-pos(k,:))+cc(2).*rand(1,D).*(gbest-pos(k,:));
    vel(k,:)=iwt.*vel(k,:)+aa(k,:); 
    vel(k,:)=(vel(k,:)>mv(k,:)).*mv(k,:)+(vel(k,:)<=mv(k,:)).*vel(k,:); 
    vel(k,:)=(vel(k,:)<(-mv(k,:))).*(-mv(k,:))+(vel(k,:)>=(-mv(k,:))).*vel(k,:);
    pos(k,:)=pos(k,:)+vel(k,:); 
    if (sum(pos(k,:)>VRmax(k,:))+sum(pos(k,:)<VRmin(k,:)))==0;
    %e(k,1)=feval(fhd,pos(k,:),func_num);
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
    
end