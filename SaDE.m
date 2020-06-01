function [DE_gbest,DE_gbestval,DE_fitcount,DE_fit_cut,DE_get_flag,ccmhist,pfithist] = SaDE(fname,...
    VTR,Max_FES,D,XRmin,XRmax,Lbound,Ubound,NP,Max_Gen,F,CR,strategy,numst,state_no,varargin)

ccmhist=[]; 
pfithist=[];
DE_get_flag = 0; 
DE_fit_cut = Max_FES;
DE_fitcount=0;

aaaa=cell(1,numst); %CR for each strategy
learngen=50;
lpcount=[];
npcount=[];
ns=[];
nf=[];
pfit=ones(1,numst);
ccm = CR*ones(1,numst);

%-----Initialize population and some arrays-------------------------------
pop = zeros(NP,D); %initialize pop to gain speed
XRRmin=repmat(XRmin,NP,1);
XRRmax=repmat(XRmax,NP,1);
rand('state',state_no);
pop=XRRmin+(XRRmax-XRRmin).*rand(NP,D);

popold    = zeros(size(pop));     % toggle population
val       = zeros(1,NP);          % create and reset the "cost array"
DE_gbest   = zeros(1,D);           % best population member ever
nfeval    = 0;                    % number of function evaluations

%------Evaluate the best member after initialization----------------------
ibest   = 1;                      % start with first population member
val(1)  = feval(fname,pop(ibest,:),varargin{:}); 
DE_gbestval = val(1);                 % best objective function value so far
nfeval  = nfeval + 1;
for i=2:NP                        % check the remaining members
    val(i) = feval(fname,pop(i,:),varargin{:}); 
    nfeval  = nfeval + 1;
    if (val(i) < DE_gbestval)           % if member is better
        ibest   = i;                 % save its location
        DE_gbestval = val(i);
    end   
end
DE_gbest = pop(ibest,:);         % best member of current iteration
ccmhist = [1,ccm];
pfithist = [1,pfit/sum(pfit)];

%------DE-Minimization---------------------------------------------
%------popold is the population which has to compete. It is--------
%------static through one iteration. pop is the newly--------------
%------emerging population.----------------------------------------

pm1 = zeros(NP,D);              % initialize population matrix 1
pm2 = zeros(NP,D);              % initialize population matrix 2
pm3 = zeros(NP,D);              % initialize population matrix 3
pm4 = zeros(NP,D);              % initialize population matrix 4
pm5 = zeros(NP,D);              % initialize population matrix 5
bm  = zeros(NP,D);              % initialize DE_gbestber  matrix
ui  = zeros(NP,D);              % intermediate population of perturbed vectors
mui = zeros(NP,D);              % mask for intermediate population
mpo = zeros(NP,D);              % mask for old population
rot = (0:1:NP-1);               % rotating index array (size NP)
rotd= (0:1:D-1);                % rotating index array (size D)
rt  = zeros(NP);                % another rotating index array
rtd = zeros(D);                 % rotating index array for exponential crossover
a1  = zeros(NP);                % index array
a2  = zeros(NP);                % index array
a3  = zeros(NP);                % index array
a4  = zeros(NP);                % index array
a5  = zeros(NP);                % index array
ind = zeros(4);


iter = 1;
while iter < Max_Gen
    popold = pop;                   % save the old population
    ind = randperm(4);              % index pointer array
    a1  = randperm(NP);             % shuffle locations of vectors
    rt = rem(rot+ind(1),NP);        % rotate indices by ind(1) positions
    a2  = a1(rt+1);                 % rotate vector locations
    rt = rem(rot+ind(2),NP);
    a3  = a2(rt+1);                
    rt = rem(rot+ind(3),NP);
    a4  = a3(rt+1);               
    rt = rem(rot+ind(4),NP);
    a5  = a4(rt+1); 
    
    pm1 = popold(a1,:);             % shuffled population 1
    pm2 = popold(a2,:);             % shuffled population 2
    pm3 = popold(a3,:);             % shuffled population 3
    pm4 = popold(a4,:);             % shuffled population 4
    pm5 = popold(a5,:);             % shuffled population 5
    
    bm = repmat(DE_gbest,NP,1);           
    
    if (iter>=learngen)
        for i=1:numst
            %             ccm(i)=median(aaaa{i}(:,1));   
            %             amax=aaaa{i}(length(aaaa{i}(:,2)),2);
            %             amin=aaaa{i}(1,2);
            %             w=(aaaa{i}(:,2)-amin)./sum(aaaa{i}(:,2)-amin);
            %             ccm(i)=sum(w.*aaaa{i}(:,1));
            %             if length(aaaa{i}(:,2))<10,
            %                 aaaa{i}(:,2)
            %             end
            %             if   ~isempty(find(aaaa{i}(:,2)==iter-1)) 
            %                 d_index=find(aaaa{i}(:,2)==aaaa{i}(1,2));
            %                 aaaa{i}(d_index,:)=[];
            %             end
            if   ~isempty(aaaa{i}) 
                ccm(i)=median(aaaa{i}(:,1));   
                d_index=find(aaaa{i}(:,2)==aaaa{i}(1,2));
                aaaa{i}(d_index,:)=[];
            else
                ccm(i)=rand;
            end
        end
    end
    
    for i=1:numst
        cc_tmp=[];
        for k=1:NP
            tt=normrnd(ccm(i),0.1);
            while tt>1 || tt<0
                tt=normrnd(ccm(i),0.1);
            end
            cc_tmp=[cc_tmp;tt];
        end
        cc(:,i)=cc_tmp;
    end
    
    % Stochastic universal sampling               %choose strategy
    rr=rand;    
    spacing=1/NP;
    randnums=sort(mod(rr:spacing:1+rr-0.5*spacing,1));  
    
    normfit=pfit/sum(pfit);
    partsum=0;   
    count(1)=0; 
    stpool=[];
    
    for i=1:length(pfit)
        partsum=partsum+normfit(i);
        count(i+1)=length(find(randnums<partsum));
        select(i,1)=count(i+1)-count(i);
        stpool=[stpool;ones(select(i,1),1)*i];
    end
    stpool = stpool(randperm(NP));
    
    for i=1:numst
        atemp=zeros(1,NP);
        aaa{i}=atemp;
        index{i}=[];
        if ~isempty(find(stpool == i))
            index{i} = find(stpool == i);
            atemp(index{i})=1;
            aaa{i}=atemp;
        end
    end
    
    aa=zeros(NP,D);
    for i=1:numst
        aa(index{i},:) = rand(length(index{i}),D) < repmat(cc(index{i},i),1,D);          % all random numbers < CR are 1, 0 otherwise
    end
    mui=aa;
    if (strategy > 1)
        st = strategy-1;		  % binomial crossover
    else
        st = strategy;		  % exponential crossover
        mui=sort(mui');	          % transpose, collect 1's in each column
        for i=1:NP
            n=floor(rand*D);
            if n > 0
                rtd = rem(rotd+n,D);
                mui(:,i) = mui(rtd+1,i); %rotate column i by n
            end
        end
        mui = mui';			  % transpose back
    end
    % jrand
    dd=ceil(D*rand(NP,1));
    for kk=1:NP
        mui(kk,dd(kk))=1; 
    end
    mpo = mui < 0.5;                % inverse mask to mui
    
%     varx=var(pop);
%     varf=var(val);
%     varfhist=[varfhist,varf];
%     for j=1:D
%         if varxold(j)<1e-12
%             c(j)=gamma;
%         else
%             c(j)=gamma*(1e+6)*varx(j)/((1e+6)*varxold(j));
%         end
%     end
    
    for i=1:numst
        %-----------jitter---------
        F=[];
        m=length(index{i});
        %         for k=1:m
        %             tt=normrnd(F0,0.3);
        %             while tt > 2 | tt < 0
        %                 tt=normrnd(F0,0.3);
        %             end
        %             F=[F;tt];
        %         end  
        %         F=repmat(F,1,D);
        %        
        F=normrnd(0.5,0.3,m,1);
%         tol=0.5;
%         F=0.8.*exp(tol.*(normrnd(0,1,m,1)-tol/2));
        F=repmat(F,1,D);
%         %         F=normrnd(0.5,0.3,m,D);    
%          if lamda==0
%             t = (1-CR).^2./NP + (NP-1)/NP;
%         else
%             t = (1-CR).^2/NP + (NP-1)/NP*(CR*(1-lamda)^2 + 1-CR + K*CR*lamda*lamda*(1-CR));
%         end
%         F = (c<t).*Fmin + (c>=t).*sqrt((c-t)./2*NP);               
%         F = (F<Fmin).*Fmin + (F>=Fmin).*F;
%         F = (F>Fmax).*Fmax + (F<=Fmax).*F;
        %        
        if i==1
            ui(index{i},:) = pm3(index{i},:) + F.*(pm1(index{i},:) - pm2(index{i},:));        % differential variation
            ui(index{i},:) = popold(index{i},:).*mpo(index{i},:) + ui(index{i},:).*mui(index{i},:);     % crossover
        end
        if i==2
            ui(index{i},:) = popold(index{i},:) + F.*(bm(index{i},:)-popold(index{i},:)) + F.*(pm1(index{i},:) - pm2(index{i},:) + pm3(index{i},:) - pm4(index{i},:));       % differential variation
            ui(index{i},:) = popold(index{i},:).*mpo(index{i},:) + ui(index{i},:).*mui(index{i},:);     % crossover
        end
        if i==3
            ui(index{i},:) = pm5(index{i},:) + F.*(pm1(index{i},:) - pm2(index{i},:) + pm3(index{i},:) - pm4(index{i},:));       % differential variation
            ui(index{i},:) = popold(index{i},:).*mpo(index{i},:) + ui(index{i},:).*mui(index{i},:);     % crossover
        end
        if i==4     
%             ----------dizzer-----------------
%                         F=[];
% %                         %             %             Fmin=sqrt((1-cc(index{i},i)/2)/NP);
% %                         %             m=length(index{i});
% %                         %             %             F=normrnd(Fmin+0.5,0.5/3);     
%                         F=repmat(normrnd(ccm(i),0.15,m,1),1,D);
%                         ui(index{i},:) = popold(index{i},:) + repmat(cc(index{i},i),1,D).*(pm5(index{i},:)-popold(index{i},:)) + repmat(cc(index{i},i),1,D).*F.*(pm1(index{i},:) - pm2(index{i},:));       % differential variation
%                         ui(index{i},:) = popold(index{i},:) + F.*(pm5(index{i},:)-popold(index{i},:)) + F.*(pm1(index{i},:) - pm2(index{i},:));       
%             ui(index{i},:) = popold(index{i},:) + repmat(rand(m,1),1,D).*(pm5(index{i},:)-popold(index{i},:)) + F.*(pm1(index{i},:) - pm2(index{i},:));       
            ui(index{i},:) = popold(index{i},:) + rand.*(pm5(index{i},:)-popold(index{i},:)) + F.*(pm1(index{i},:) - pm2(index{i},:));  
        end
    end    
    
    for i=1:NP
        outbind=find(ui(i,:) < Lbound);
        if size(outbind,2)~=0
            %                 % Periodica
            %                 ui(i,outbind)=2*XRmin(outbind)-ui(i,outbind);
            % Random
            ui(i,outbind)=XRmin(outbind)+(XRmax(outbind)-XRmin(outbind)).*rand(1,size(outbind,2));
            %                 % Fixed
            %                 ui(i,outbind)=XRmin(outbind);
        end            
        outbind=find(ui(i,:) > Ubound);
        if size(outbind,2)~=0
            %                 % Periodica
            %                 ui(i,outbind)=2*XRmax(outbind)-ui(i,outbind);
            % Random
            ui(i,outbind)=XRmin(outbind)+(XRmax(outbind)-XRmin(outbind)).*rand(1,size(outbind,2));
            %                 % Fixed
            %                 ui(i,outbind)=XRmax(outbind);
        end
    end
    lpcount=zeros(1,numst); 
    npcount=zeros(1,numst);
    for i=1:NP
        tempval = feval(fname,ui(i,:),varargin{:});   % check cost of competitor
        nfeval  = nfeval + 1;
        if (tempval <= val(i))  % if competitor is better than value in "cost array"
            pop(i,:) = ui(i,:);  % replace old vector with new one (for new iteration)
            val(i)   = tempval;  % save value in "cost array"
            tlpcount=zeros(1,numst);
            for j=1:numst
                temp=aaa{j};
                tlpcount(j)=temp(i);
                if tlpcount(j)==1
                    aaaa{j}=[aaaa{j};cc(i,j) iter];  
                end
            end
            lpcount=[lpcount;tlpcount];
            %----we update DE_gbestval only in case of success to save time-----------
            if (tempval <= DE_gbestval)     % if competitor better than the best one ever
                DE_gbestval = tempval;      % new best value
                DE_gbest = ui(i,:);      % new best parameter vector ever
                if DE_gbestval <= VTR & DE_get_flag == 0
                    DE_fit_cut=nfeval;
                    DE_get_flag=1;
%                     DE_fitcount = nfeval;
%                     return;
                end
            end
        else
            tnpcount=zeros(1,numst);
            for j=1:numst
                temp=aaa{j};
                tnpcount(j)=temp(i);
            end
            npcount=[npcount;tnpcount];
        end
        
        if nfeval+1 > Max_FES
            DE_fitcount = Max_FES;
            pfithist = [pfithist;[iter+2,pfit/sum(pfit)]];
            ccmhist = [ccmhist;[iter+2,ccm]];
            return;
        end
        %     end
    end %---end for imember=1:NP
    pfithist = [pfithist;[iter+2,pfit/sum(pfit)]];
    ccmhist = [ccmhist;[iter+2,ccm]];    
    ns=[ns;sum(lpcount,1)];
    nf=[nf;sum(npcount,1)];    
    
    if iter >= learngen,
        %         ww=repmat((1:learngen)',1,numst);
        %         ns=ww.*ns;
        %         nf=ww.*nf;
        for i=1:numst            
            if (sum(ns(:,i))+sum(nf(:,i))) == 0
                pfit(i) = 0.01;
            else
                pfit(i) = sum(ns(:,i))/(sum(ns(:,i))+ sum(nf(:,i))) + 0.01;
            end
        end
        if ~isempty(ns), ns(1,:)=[];   end
        if ~isempty(nf), nf(1,:)=[];   end
    end 
    iter = iter + 1;
%     if rem(iter,50)==0
%         nfeval
%         DE_gbestval
%         %         DE_gbest

% %         plot(pop(:,1),pop(:,2),'r.');
% %         hold on
% %         plot(DE_gbest(1),DE_gbest(2),'bo')
% %         hold off
% %         drawnow
%         %     %         sum(ns(2,:))+sum(nf(2,:))
%         %             size(ns)
%         %             size(nf)
%     end
end %---end while ((iter < Max_Gen) ...