function [zyellow, zred, both, no] = generate(p, num_sim, pop_size, init_copy_yellow, init_copy_red, s)

Zred=[];
Zyellow=[];
Both=[];
No=[];


%%PLOT OF BINOMIAL RANDOM NUMBERS WITH AND WITHOUT SELECTION
for N=1:num_sim
    

Ryellow=[init_copy_yellow];
Rred=[init_copy_red];

if init_copy_red>0 && init_copy_yellow>0
    redplot=[0];
    yellowplot=[0];
    mixplot=[1];
    noplot=[0];
    sumy=0;
    sumr=0;
    summ=1;
    sumno=0;
elseif init_copy_yellow>0 && init_copy_red==0
    redplot=[0];
    yellowplot=[1];
    mixplot=[0];
    noplot=[0];
    sumy=1;
    sumr=0;
    summ=0;
    sumno=0;
elseif init_copy_yellow==00 && init_copy_red>0
    redplot=[1];
    yellowplot=[0];
    mixplot=[0];
    noplot=[0];
    sumy=0;
    sumr=1;
    summ=0;
    sumno=0;
else
    redplot=[0];
    yellowplot=[0];
    mixplot=[0];
    noplot=[1];
    sumy=0;
    sumr=0;
    summ=0;
    sumno=1;
end
    


while length(Ryellow)<pop_size
  
    idx=find(Ryellow==0 & Rred==0);
    no=length(idx); %%count cells with any copy
    yes=length(Ryellow)-no; %%count cells with at least 1 ecDNA copy
    xi1=rand;
    xi2=rand;
    %%tauy=1/(2*yellow)*log(1/xi1);
    %%taur=1/(red)*log(1/xi1);
    tauy=random('Exponential',1/(s*yes));
    taur=random('Exponential',1/(1*no));


    if tauy<taur % ecDNA positive divides
        j=random("Discrete uniform", length(Ryellow));
        while Ryellow(j)==0 & Rred(j)==0
          j=random("Discrete uniform", length(Ryellow));
        end
        
        if Ryellow(j)~=0 & Rred(j)~=0
            summ=summ-1;
        elseif Ryellow(j)==0 & Rred(j)~=0
            sumr=sumr-1;
        elseif Ryellow(j)~=0 & Rred(j)==0
            sumy=sumy-1;
        else
            sumno=sumno-1;
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
    %first daughter
    Ryellow=[Ryellow,appob1+appob3];
    Rred=[Rred,appor1+appor3];
    %second daughter
    Ryellow=[Ryellow,appob2+appob4];
    Rred=[Rred,appor2+appor4];
    

    Ryellow(j)=[]; %mother is eliminated
    Rred(j)=[];
    
    %adding daughters to trends
    if Ryellow(end)==0 && Rred(end)~=0
        sumr=sumr+1;
    elseif Ryellow(end)~=0 && Rred(end)==0
        sumy=sumy+1;
    elseif Ryellow(end)~=0 && Rred(end)~=0
        summ=summ+1;
    else
        sumno=sumno+1;
    end

    if Ryellow(end-1)==0 && Rred(end-1)~=0
        sumr=sumr+1;
    elseif Ryellow(end-1)~=0 && Rred(end-1)==0
        sumy=sumy+1;
    elseif Ryellow(end-1)~=0 && Rred(end-1)~=0
        summ=summ+1;
    else
        sumno=sumno+1;
    end
    
    
    else  %ecDNA negative divides
    sumno=sumno+1;
    Ryellow=[Ryellow,0];
    Rred=[Rred,0];

    end

    
    
    redplot=[redplot,sumr/length(Ryellow)];
    yellowplot=[yellowplot,sumy/length(Ryellow)];
    mixplot=[mixplot,summ/length(Ryellow)];
    noplot=[noplot,sumno/length(Ryellow)];
  
end

Zyellow=[Zyellow;redplot];
Zred=[Zred;yellowplot];
Both=[Both;mixplot];
No=[No;noplot];
   

end

zred=mean(Zred);
zyellow=mean(Zyellow);
both=mean(Both);
no=mean(No);


end
