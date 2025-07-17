%% Fraction of subpopulation over time for different starting points of switching / EXAMPLE WITH NEUTRAL SELECTION

clear all

%Parameters
p=0.01;  %p_y=p_r
num_sim=1000;
pop_size=100000;
init_copy_yellow=1;
init_copy_red=0;


%SWITCHING AT 0

[yellow_0, red_0, both_0] = generate_later_switching_neutral(p, num_sim, pop_size, init_copy_yellow, init_copy_red, 0);


%SWITCHING AT 10
[yellow_10, red_10, both_10] = generate_later_switching_neutral(p, num_sim, pop_size, init_copy_yellow, init_copy_red, 10);


%SWITCHING AT 100
[yellow_100, red_100, both_100] = generate_later_switching_neutral(p, num_sim, pop_size, init_copy_yellow, init_copy_red, 100);


f1=figure;

y=linspace(1,length(yellow_0),length(yellow_0));

semilogx(y,yellow_0,'y');
hold on
semilogx(y,red_0,'r');
hold on
semilogx(y,both_0,'m');
hold on
semilogx(y,yellow_10,'y--');
hold on
semilogx(y,red_10,'r--');
hold on
semilogx(y,both_10,'m--');
hold on
semilogx(y,yellow_100,'y-.');
hold on
semilogx(y,red_100,'r-.');
hold on
semilogx(y,both_100,'m-.');







xlabel("Number of cells");
ylabel("Frequency");


function [zyellow, zred, both] = generate_later_switching_neutral(p, num_sim, pop_size, init_copy_yellow, init_copy_red, t_s)

Zred=[];
Zyellow=[];
Both=[];


for N=1:num_sim
    

Ryellow=[init_copy_yellow];
Rred=[init_copy_red];
if init_copy_red>0 && init_copy_yellow>0
    sumy=0;
    sumr=0;
    summ=1;
    redplot=[0];
    yellowplot=[0];
    mixplot=[1];
elseif init_copy_yellow>0 && init_copy_red==0
    sumy=1;
    sumr=0;
    summ=0;
    redplot=[0];
    yellowplot=[1];
    mixplot=[0];
elseif init_copy_yellow==0 && init_copy_red>0
    sumy=0;
    sumr=1;
    summ=0;
    redplot=[1];
    yellowplot=[0];
    mixplot=[0];
end


if t_s>0

    while length(Ryellow) < t_s
    
        j=random("Discrete uniform", length(Ryellow));  
        if Ryellow(j)~=0 & Rred(j)~=0
           summ=summ-1;
        elseif Ryellow(j)==0 & Rred(j)~=0
            sumr=sumr-1;
        elseif Ryellow(j)~=0 & Rred(j)==0
            sumy=sumy-1;
        end
    
        copyb=binornd(2*Ryellow(j),1/2);
        Ryellow=[Ryellow,copyb,2*Ryellow(j)-copyb];
        Rred=[Rred,0,0];
    
        Ryellow(j)=[];
        Rred(j)=[];
        
        if Ryellow(end)==0 && Rred(end)~=0
            sumr=sumr+1;
        elseif Ryellow(end)~=0 && Rred(end)==0
            sumy=sumy+1;
        elseif Ryellow(end)~=0 && Rred(end)~=0
            summ=summ+1;
        end
    
        if Ryellow(end-1)==0 && Rred(end-1)~=0
            sumr=sumr+1;
        elseif Ryellow(end-1)~=0 && Rred(end-1)==0
            sumy=sumy+1;
        elseif Ryellow(end-1)~=0 && Rred(end-1)~=0
            summ=summ+1;
        end
    
        redplot=[redplot,sumr/length(Ryellow)];
        yellowplot=[yellowplot,sumy/length(Ryellow)];
        mixplot=[mixplot,summ/length(Ryellow)];
    
    
    
        %noswitching
    end 
end

while length(Ryellow)<pop_size
  
    j=random("Discrete uniform", length(Ryellow));
    %removal of mother cell from trends
    if Ryellow(j)~=0 & Rred(j)~=0
       summ=summ-1;
    elseif Ryellow(j)==0 & Rred(j)~=0
        sumr=sumr-1;
    elseif Ryellow(j)~=0 & Rred(j)==0
        sumy=sumy-1;
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
    %Rred=[Rred,copyr+appor];   %%aggiungo copie rosse
    
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
    %Rred=[Rred,2*Rred(j)-copyr+appor];   %%aggiungo copie rosse
    
    
    appob3=0;
    appor3=copyr;
    if copyr~=0
      for h=1:copyr
        xi=rand;
        if xi<p   %PROBABILITA ROSSA
            appob3=appob3+1;
            appor3=appor3-1;
        end
    end
    end
    %Ryellow(end-1)=appob;
    %Rred(end-1)=appor;   %%aggiungo copie rosse
    
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
    

    Ryellow(j)=[];  %removal mother
    Rred(j)=[];
    
    %adding first daughter to pure yellow trend if yellow, pure red trend
    %if red or mixed if mixed 
    if Ryellow(end)==0 && Rred(end)~=0
        sumr=sumr+1;
    elseif Ryellow(end)~=0 && Rred(end)==0
        sumy=sumy+1;
    elseif Ryellow(end)~=0 && Rred(end)~=0
        summ=summ+1;
    end

    if Ryellow(end-1)==0 && Rred(end-1)~=0
        sumr=sumr+1;
    elseif Ryellow(end-1)~=0 && Rred(end-1)==0
        sumy=sumy+1;
    elseif Ryellow(end-1)~=0 && Rred(end-1)~=0
        summ=summ+1;
    end

    redplot=[redplot,sumr/length(Ryellow)];
    yellowplot=[yellowplot,sumy/length(Ryellow)];
    mixplot=[mixplot,summ/length(Ryellow)];

end

Zyellow=[Zyellow;yellowplot];
Zred=[Zred;redplot];
Both=[Both;mixplot];
   

end

zred=mean(Zred);
zyellow=mean(Zyellow);
both=mean(Both);


end

