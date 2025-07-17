%% Distribution of ecDNA copies among different subpopulations

clear all

%Parameters
p=0.99;  %p_y=p_r
num_sim=1000;
pop_size=100000;
init_copy_yellow=1;
init_copy_red=0;

%%nonid=0 means identical selection, nonid=1 means non identical selection


%NEUTRAL SELECTION CASE

[yellowblockneu, redblockneu, mixblockneu]=generate_distribution(p, num_sim, pop_size, init_copy_yellow, init_copy_red, 1, 0);

%IDENTICAL SELECTION

[yellowblocksel, redblocksel, mixblocksel]=generate_distribution(p, num_sim, pop_size, init_copy_yellow, init_copy_red, 2, 0);

%NON IDENTICAL SELECTION

[yellowblockselnonid, redblockselnonid, mixblockselnonid]=generate_distribution(p, num_sim, pop_size, init_copy_yellow, init_copy_red, 2, 1);



%save('simulation_distributionpuremix099.mat');




%Plotting commands


xlabel("Copies of ecDNA");
ylabel("Frequency");

x1=linspace(1,max(yellowblockneu),max(yellowblockneu));
f1 = ksdensity(yellowblockneu,x1);
sem1 = std(yellowblockneu) / sqrt(length(yellowblockneu)); 
x2=linspace(1,max(redblockneu),max(redblockneu));
f2 = ksdensity(redblockneu,x2);
sem2 = std(redblockneu) / sqrt(length(redblockneu)); 
x3=linspace(1,max(mixblockneu),max(mixblockneu));
f3 = ksdensity(mixblockneu,x3);
sem3 = std(mixblockneu) / sqrt(length(mixblockneu)); 

%Density curves
hold on
plot(x1, f1, 'y', 'LineWidth', 1.5); 
hold on
plot(x2, f2, 'r', 'LineWidth', 1.5); 
hold on
plot(x3, f3, 'm', 'LineWidth', 1.5); 

%Error bars (opt)
%hold on
%errorbar(x1, f1, sem1 * ones(size(f1)), 'y.', 'MarkerSize', 10); 
%hold on
%errorbar(x2, f2, sem2 * ones(size(f2)), 'r.', 'MarkerSize', 10); 
%hold on
%errorbar(x3, f3, sem3 * ones(size(f3)), 'm.', 'MarkerSize', 10); 


x1=linspace(1,max(yellowblocksel),max(yellowblocksel));
f1 = ksdensity(yellowblocksel,x1);
sem1 = std(yellowblocksel) / sqrt(length(yellowblocksel)); 
x2=linspace(1,max(redblocksel),max(redblocksel));
f2 = ksdensity(redblocksel,x2);
sem2 = std(redblocksel) / sqrt(length(redblocksel)); 
x3=linspace(1,max(mixblocksel),max(mixblocksel));
f3 = ksdensity(mixblocksel,x3);
sem3 = std(mixblocksel) / sqrt(length(mixblocksel)); 

%Density curves
hold on
plot(x1, f1, 'y--', 'LineWidth', 1.5); 
hold on
plot(x2, f2, 'r--', 'LineWidth', 1.5); 
hold on
plot(x3, f3, 'm--', 'LineWidth', 1.5);

%Error bars (opt)
%hold on
%errorbar(x1, f1, sem1 * ones(size(f1)), 'y.', 'MarkerSize', 10);
%hold on
%errorbar(x2, f2, sem2 * ones(size(f2)), 'r.', 'MarkerSize', 10); 
%hold on
%errorbar(x3, f3, sem3 * ones(size(f3)), 'm.', 'MarkerSize', 10); 

x1=linspace(1,max(yellowblockselnonid),max(yellowblockselnonid));
f1 = ksdensity(yellowblockselnonid,x1);
sem1 = std(yellowblockselnonid) / sqrt(length(yellowblockselnonid)); 
x2=linspace(1,max(redblockselnonid),max(redblockselnonid));
f2 = ksdensity(redblockselnonid,x2);
sem2 = std(redblockselnonid) / sqrt(length(redblockselnonid)); 
x3=linspace(1,max(mixblockselnonid),max(mixblockselnonid));
f3 = ksdensity(mixblockselnonid,x3);
sem3 = std(mixblockselnonid) / sqrt(length(mixblockselnonid)); 

%Density curves
hold on
plot(x1, f1, 'y-.', 'LineWidth', 1.5);
hold on
plot(x2, f2, 'r-.', 'LineWidth', 1.5); 
hold on
plot(x3, f3, 'm-.', 'LineWidth', 1.5); 

%Error bars (opt)
%hold on
%errorbar(x1, f1, sem1 * ones(size(f1)), 'y.', 'MarkerSize', 10); 
%hold on
%errorbar(x2, f2, sem2 * ones(size(f2)), 'r.', 'MarkerSize', 10); 
%hold on
%errorbar(x3, f3, sem3 * ones(size(f3)), 'm.', 'MarkerSize', 10); 








function [yellowblock, redblock, mixblock] = generate_distribution(p, num_sim, pop_size, init_copy_yellow, init_copy_red, s, nonid)
Ared_n=[]; %neutral
Ayellow_n=[];
Ared_s=[]; %identical selection
Ayellow_s=[];
Ared_ns=[]; %non identical selection
Ayellow_ns=[];

if s==1  %NEUTRAL SELECTION CASE
    
    for N=1:num_sim

    Ryellow=[init_copy_yellow];
    Rred=[init_copy_red];
    
    while length(Ryellow)<pop_size
      
        j=random("Discrete uniform", length(Ryellow));
        copyr=binornd(2*Rred(j),1/2);
        copyb=binornd(2*Ryellow(j),1/2);   %binomial allocation
        
        %%switching
        
        appob1=copyb;
        appor1=0;
        if copyb~=0
          for h=1:copyb
            xi=rand;
            if xi<p
                appob1=appob1-1;
                appor1=appor1+1;
            end
        end
        end
        %Ryellow=[Ryellow,appob];
        %Rred=[Rred,copyr+appor];   
        
        appob2=2*Ryellow(j)-copyb;
        appor2=0;
        if copyb~=2*Ryellow(j)
           for h=1:(2*Ryellow(j)-copyb)
            xi=rand;
            if xi<p
                appob2=appob2-1;
                appor2=appor2+1;
            end
        end
        end
        %Ryellow=[Ryellow,appob];
        %Rred=[Rred,2*Rred(j)-copyr+appor];   
        
        
        appob3=0;
        appor3=copyr;
        if copyr~=0
          for h=1:copyr
            xi=rand;
            if xi<p   
                appob3=appob3+1;
                appor3=appor3-1;
            end
        end
        end
        %Ryellow(end-1)=appob;
        %Rred(end-1)=appor;   
        
        appor4=2*Rred(j)-copyr;
        appob4=0;
        if 2*Rred(j)-copyr~=0
           for h=1:2*Rred(j)-copyr
            xi=rand;
            if xi<p
                appob4=appob4+1;
                appor4=appor4-1;
            end
        end
        end
        
        %first daughter
        Ryellow=[Ryellow,appob1+appob3]; 
        Rred=[Rred,appor1+appor3];
        %second daughter
        Ryellow=[Ryellow,appob2+appob4];
        Rred=[Rred,appor2+appor4];
    
    
        Ryellow(j)=[];  %mother cell eliminated
        Rred(j)=[];
    
    
    
    end    
       
    
    
    Ared_n=[Ared_n,Rred];
    Ayellow_n=[Ayellow_n,Ryellow];
    
    end
    
    
    idxb=find(Ared_n==0);  %pure yellow
    idxr=find(Ayellow_n==0);  %pure red
    yellowblock=Ayellow_n(idxb);   
    redblock=Ared_n(idxr);   
    idxm=find(Ared_n~=0 & Ayellow_n~=0); %mixed
    mixblock=Ared_n(idxm)+Ayellow_n(idxm);

elseif s>1 && nonid==0  %IDENTICAL SELECTION CASE
    
    for N=1:num_sim

    Ryellow=[init_copy_yellow];
    Rred=[init_copy_red];
    while length(Ryellow)<pop_size
        idx=find(Ryellow==0 & Rred==0);
        no=length(idx); %%count cells no copies
        yes=length(Ryellow)-no; %%count cells with at least 1 copy of any colour
        xi1=rand;
        xi2=rand;
        %%tauy=1/(2*yellow)*log(1/xi1);
        %%taur=1/(red)*log(1/xi1);
        tauyes=random('Exponential',1/(s*yes));
        tauno=random('Exponential',1/(1*no));  %Gillespie
    
    
        if tauyes<tauno    %ecDNA positive divides
            j=random("Discrete uniform", length(Ryellow));
            while Ryellow(j)==0 && Rred(j)==0
              j=random("Discrete uniform", length(Ryellow));
            end
        copyr=binornd(2*Rred(j),1/2);
        copyb=binornd(2*Ryellow(j),1/2);   %binomial allocation
        
        appob1=copyb;
        appor1=0;
        if copyb~=0
          for h=1:copyb
            xi=rand;
            if xi<p
                appob1=appob1-1;
                appor1=appor1+1;
            end
        end
        end
        %Ryellow=[Ryellow,appob];
        %Rred=[Rred,copyr+appor];  
        
        appob2=2*Ryellow(j)-copyb;
        appor2=0;
        if copyb~=2*Ryellow(j)
           for h=1:(2*Ryellow(j)-copyb)
            xi=rand;
            if xi<p
                appob2=appob2-1;
                appor2=appor2+1;
            end
        end
        end
        %Ryellow=[Ryellow,appob];
        %Rred=[Rred,2*Rred(j)-copyr+appor];   
        
        
        appob3=0;
        appor3=copyr;
        if copyr~=0
          for h=1:copyr
            xi=rand;
            if xi<p   
                appob3=appob3+1;
                appor3=appor3-1;
            end
        end
        end
        %Ryellow(end-1)=appob;
        %Rred(end-1)=appor;  
        
        appor4=2*Rred(j)-copyr;
        appob4=0;
        if 2*Rred(j)-copyr~=0
           for h=1:2*Rred(j)-copyr
            xi=rand;
            if xi<p
                appob4=appob4+1;
                appor4=appor4-1;
            end
        end
        end
    
        %first daughter
        Ryellow=[Ryellow,appob1+appob3];
        Rred=[Rred,appor1+appor3];
        %second daughter
        Ryellow=[Ryellow,appob2+appob4];
        Rred=[Rred,appor2+appor4];
    
    
        Ryellow(j)=[]; %mother cell eliminated
        Rred(j)=[];
        
        else  %ecDNA negative divides
    
        Ryellow=[Ryellow,0];
        Rred=[Rred,0];
    
        end
    
    
    end    
       
    
    
    Ared_s=[Ared_s,Rred];
    Ayellow_s=[Ayellow_s,Ryellow];
    
    end
    
    
    idxb=find(Ared_s==0);  %pure yellow
    idxr=find(Ayellow_s==0);  %pure red 
    yellowblock=Ayellow_s(idxb);   
    redblock=Ared_s(idxr);   
    idxm=find(Ared_s~=0 & Ayellow_s~=0); %mixed
    mixblock=Ared_s(idxm)+Ayellow_s(idxm);
    %idxn=find(Ared==0 & Ayellow==0);   

elseif s>1 && nonid==1  %NON IDENTICAL SELECTION CASE
    
    for N=1:num_sim

    Ryellow=[init_copy_yellow];
    Rred=[init_copy_red];
    
    while length(Ryellow)<pop_size
        idx=find(Ryellow~=0);
        yes=length(idx); %%count cells with either red copies or no copies
        no=length(Ryellow)-yes; %%count cells with at least 1 yellow copy
        xi1=rand;
        xi2=rand;
        %%tauy=1/(2*yellow)*log(1/xi1);
        %%taur=1/(red)*log(1/xi1);
        tauy=random('Exponential',1/(s*yes));
        taur=random('Exponential',1/(1*no));
    
    
        if tauy<taur  %either pure yellow or mixed divide
           j=random("Discrete uniform", length(Ryellow));
            while Ryellow(j)==0 && Rred(j)==0
              j=random("Discrete uniform", length(Ryellow));
            end
        copyr=binornd(2*Rred(j),1/2);
        copyb=binornd(2*Ryellow(j),1/2);
        
        appob1=copyb;
        appor1=0;
        if copyb~=0
          for h=1:copyb
            xi=rand;
            if xi<p
                appob1=appob1-1;
                appor1=appor1+1;
            end
        end
        end
        %Ryellow=[Ryellow,appob];
        %Rred=[Rred,copyr+appor];  
        
        appob2=2*Ryellow(j)-copyb;
        appor2=0;
        if copyb~=2*Ryellow(j)
           for h=1:(2*Ryellow(j)-copyb)
            xi=rand;
            if xi<p
                appob2=appob2-1;
                appor2=appor2+1;
            end
        end
        end
        %Ryellow=[Ryellow,appob];
        %Rred=[Rred,2*Rred(j)-copyr+appor];   
        
        
        appob3=0;
        appor3=copyr;
        if copyr~=0
          for h=1:copyr
            xi=rand;
            if xi<p  
                appob3=appob3+1;
                appor3=appor3-1;
            end
        end
        end
        %Ryellow(end-1)=appob;
        %Rred(end-1)=appor;   
        
        appor4=2*Rred(j)-copyr;
        appob4=0;
        if 2*Rred(j)-copyr~=0
           for h=1:2*Rred(j)-copyr
            xi=rand;
            if xi<p
                appob4=appob4+1;
                appor4=appor4-1;
            end
        end
        end
        Ryellow=[Ryellow,appob1+appob3];
        Rred=[Rred,appor1+appor3];
        Ryellow=[Ryellow,appob2+appob4];
        Rred=[Rred,appor2+appor4];
    
    
        Ryellow(j)=[];
        Rred(j)=[];
        
        else  %either pure red or ecDNA negative divide
            j=random("Discrete uniform", length(Ryellow));
            while Ryellow(j)~=0
              j=random("Discrete uniform", length(Ryellow));
            end
        copyr=binornd(2*Rred(j),1/2); 
        appob=0;
        appor=copyr;
        if copyr~=0
          for h=1:copyr
            xi=rand;
            if xi<p   
                appob=appob+1;
                appor=appor-1;
            end
        end
        end
        Ryellow=[Ryellow,appob];
        Rred=[Rred,appor];   
        
        appor=2*Rred(j)-copyr;
        appob=0;
        if copyr~=2*Rred(j)
           for h=1:(2*Rred(j)-copyr)
            xi=rand;
            if xi<p
                appob=appob+1;
                appor=appor-1;
            end
        end
        end
        Ryellow=[Ryellow,appob];
        Rred=[Rred,appor];
    
        Ryellow(j)=[];
        Rred(j)=[];
    
        end
    
    
    end    
       
    
    
    Ared_ns=[Ared_ns,Rred];
    Ayellow_ns=[Ayellow_ns,Ryellow];
    
    end
    
    
    idxb=find(Ared_ns==0);  %pure yellow
    idxr=find(Ayellow_ns==0);  %pure red
    yellowblock=Ayellow_ns(idxb);   
    redblock=Ared_ns(idxr);   
    idxm=find(Ared_ns~=0 & Ayellow_ns~=0); %mixed
    mixblock=Ared_ns(idxm)+Ayellow_ns(idxm);
end

end

















