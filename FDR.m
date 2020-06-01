function [position,value,iteration]= FDR(benchmark_func_org,Particle_Number,Dimension,Max_Gen,Max_FES,func_num,runs,range)
%[gbest,gbestval,fitcount]= PSO_func('f8',3500,200000,30,30,-5.12,5.12)

% range=[-100,100];

% result_data=zeros(1,Max_FES);
% result=zeros(1,Particle_Number); % Initial fitness values
% check_vel=[];
me=Max_Gen;
ps=Particle_Number;
D=Dimension;
outcome=[];
Pg=zeros(1,Max_FES); 
fii=[1 1 2];
y=(-1400:100:1400);
% load fbias_data
% y=f_bias;
iwt=0.9-(1:me).*(0.5./me);
for time=1:runs
    rand('state',sum(100*clock));
    fitold=1;
    VRmin=range(1);
    VRmax=range(2);
    if length(VRmin)==1
        VRmin=repmat(VRmin,1,D);
        VRmax=repmat(VRmax,1,D);
    end
    mv=0.2*(VRmax-VRmin);
    VRmin=repmat(VRmin,ps,1);
    VRmax=repmat(VRmax,ps,1);
    Vmin=repmat(-mv,ps,1);
    Vmax=-Vmin;
    pos=VRmin+(VRmax-VRmin).*rand(ps,D);
    fitcount=0;
    % for i=1:ps;
         result=benchmark_func_org(pos',func_num);
         fitcount=fitcount+ps;
    %      result_data(fitcount)=result(i);
    % end

    vel=Vmin+2.*Vmax.*rand(ps,D);%initialize the velocity of the particles
    pbest=pos;
    pbestval=result; %initialize the pbest and the pbest's fitness value
    [gbestval,gbestid]=min(pbestval);
    gbest=pbest(gbestid,:);%initialize the gbest and the gbest's fitness value
    gbestrep=repmat(gbest,ps,1);
    i=0;

    while i<me
          i=i+1;

        for k=1:ps
            dis=abs(repmat(pbest(k,:),ps,1)-pbest(1:ps,:));
            fiterr=repmat(pbestval(k),ps,1)-pbestval(1:ps)';
            fiterr=repmat(fiterr,1,D);
            fiterr=fiterr-(dis==zeros(ps,D)).*fiterr;
            dis=dis+(dis==zeros(ps,D));
            FDR=fiterr./dis;
            [fdr,Fid]=max(FDR);
            for dimcnt=1:D
                Pnd(k,dimcnt)=pbest(Fid(dimcnt),dimcnt);
            end
            aa(k,:)=fii(1).*rand(1,D).*(pbest(k,:)-pos(k,:))+fii(2).*rand(1,D).*(gbestrep(k,:)-pos(k,:))+fii(3).*rand(1,D).*(Pnd(k,:)-pos(k,:));
            vel(k,:)=iwt(i).*vel(k,:)+aa(k,:); 
            vel(k,:)=(vel(k,:)>mv).*mv+(vel(k,:)<=mv).*vel(k,:); 
            vel(k,:)=(vel(k,:)<(-mv)).*(-mv)+(vel(k,:)>=(-mv)).*vel(k,:);
            pos(k,:)=pos(k,:)+vel(k,:); 

            if (sum(pos(k,:)>VRmax(k,:))+sum(pos(k,:)<VRmin(k,:)))==0;
                result(k)=benchmark_func_org(pos(k,:)',func_num);
                fitcount=fitcount+1;
    %             result_data(fitcount)=result(k);
                if fitcount>=Max_FES
                    break;
                end
                tmp=(pbestval(k)<result(k));
                temp=repmat(tmp,1,D);
                pbest(k,:)=temp.*pbest(k,:)+(1-temp).*pos(k,:);
                pbestval(k)=tmp.*pbestval(k)+(1-tmp).*result(k);%update the pbest
                if pbestval(k)<gbestval
                   gbest=pbest(k,:);
                   gbestval=pbestval(k);
                   gbestrep=repmat(gbest,ps,1);%update the gbest
                end
           end
        end
        if fitcount>=Max_FES
           break;
        end
        if time==1
            Pg(1,fitold:fitcount)=gbestval-y(func_num); 
        end
        fitold=fitcount;
    if (i==me)&&(fitcount<Max_FES)
        i=i-1;
    end
   end
    outcome=[outcome gbestval-y(func_num)];
end
position=Pg;
value=outcome;
iteration=k;
end


