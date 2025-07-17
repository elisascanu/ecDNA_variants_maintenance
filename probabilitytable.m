%% Probability of staying in different states
    % it creates a table of probabilities
    % Column one counts mother cells that are mix
    % Column two counts mother cells that are pure
    % Column three counts mother cells that are ecDNA negative
    % Columns four to six count daughter cells from mix mothers, based on the nature of them
    % Columns seven to nine count daughter cells from pure mothers, based on the nature of them
    % Column ten records daughters of ecDNA negative mothers
    % This table can be easily used to produce the probability plots in Figure 5


clear all

%Parameters
pop_size=100000;
num_sim=1000;
p=0.8;
s=2;
init_copy_yellow=1;
init_copy_red=0;

%Initialising
C=zeros(pop_size,10);

for N=1:num_sim
    
comp1=0;
comp2=0;
comp3=0;
comp4=0;
comp5=0;
comp6=0;
Ryellow=[init_copy_yellow]; %%vector for allocating # of red ecDNA copies in each simulation
Rred=[init_copy_red]; %%vector for allocating # of orange ecDNA copies in each simulation

if init_copy_red>0 && init_copy_yellow>0
    mix=1;
    pure=0;
    nocop=0;
elseif init_copy_yellow==0 && init_copy_red==0
    mix=0;
    pure=0;
    nocop=1;
else 
    mix=0;
    pure=1;
    nocop=0;
end

somma=0;
nocop=[nocop];
while length(Ryellow)<pop_size
    idx=find(Ryellow==0 & Rred==0);
    no=length(idx); %%count cells with no ecDNA copies
    yes=length(Ryellow)-no; %%count cells with at least 1 ecDNA copy
    xi1=rand;
    xi2=rand;
    %%tauy=1/(2*yellow)*log(1/xi1);
    %%taur=1/(red)*log(1/xi1);
    tauy=random('Exponential',1/(s*yes));
    taur=random('Exponential',1/(1*no));


    if tauy<taur  %ecDNA positive divides
        j=random("Discrete uniform", length(Ryellow));
        while Ryellow(j)==0 && Rred(j)==0
          j=random("Discrete uniform", length(Ryellow));
        end
    copyr=binornd(2*Rred(j),1/2);
    copyb=binornd(2*Ryellow(j),1/2);
    
    
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
    Ryellow=[Ryellow,appob1+appob3];
    Rred=[Rred,appor1+appor3];
    Ryellow=[Ryellow,appob2+appob4];
    Rred=[Rred,appor2+appor4];
    
    appomix=0;
    appopure=0;
    appono=0;

    %aidentify nature of daughters
       if Rred(end-1)~=0 & Ryellow(end-1)~=0
           appomix=appomix+1;
       elseif Rred(end-1)==0 & Ryellow(end-1)==0
           appono=appono+1;
       else
           appopure=appopure+1;
       end

       if Rred(end)~=0 & Ryellow(end)~=0
           appomix=appomix+1;
       elseif Rred(end)==0 & Ryellow(end)==0
           appono=appono+1;
       else
           appopure=appopure+1;
       end
    
       %adding mother and daughters to the trends
        if Rred(j)~=0 && Ryellow(j)~=0
            C(length(Ryellow)-2,1)=C(length(Ryellow)-2,1)+1;
            C(length(Ryellow)-2,4)=C(length(Ryellow)-2,4)+appomix;
            C(length(Ryellow)-2,5)=C(length(Ryellow)-2,5)+appopure;
            C(length(Ryellow)-2,6)=C(length(Ryellow)-2,6)+appono;
        else
            C(length(Ryellow)-2,2)=C(length(Ryellow)-2,2)+1;
            C(length(Ryellow)-2,7)=C(length(Ryellow)-2,7)+appomix;
            C(length(Ryellow)-2,8)=C(length(Ryellow)-2,8)+appopure;
            C(length(Ryellow)-2,9)=C(length(Ryellow)-2,9)+appono;
            
        end
    
    
    
    Ryellow(j)=[]; %mother removed
    Rred(j)=[];

    else  %ecDNA negative divides
        Ryellow=[Ryellow,0];
        Rred=[Rred,0];
        C(length(Ryellow)-1,3)=C(length(Ryellow)-1,3)+1;
        C(length(Ryellow)-1,10)=C(length(Ryellow)-1,10)+2;
    end


end
end

T = table(C(:,1),C(:,2),C(:,3),C(:,4),C(:,5),C(:,6),C(:,7),C(:,8),C(:,9),C(:,10),'VariableNames',["Mother mix", "Mother pure","Mother noecDNA", "Daughter mix from mix", "Daughter pure from mix", "Daughter noecDNA from mix", "Daughter mix from pure", "Daughter pure from pure", "Daughter noecDNA from pure", "Daughter noecDNA from noecDNA"]);

writetable(T,sprintf('table_output.csv'));