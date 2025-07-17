%% Shannon index at single cell level

clear all


%Parameters
num_sim=1000;
pop_size=10000;
prob=0.01;
s=1;
init_copy_yellow=1;
init_copy_red=1;


%Initialise
H=[];
K=[];
K1=[];
Hred=[];
Kred=[];


for N=1:num_sim
    

Ryellow=[init_copy_yellow]; %%vector for allocating # of red ecDNA copies in each simulation
Rred=[init_copy_red]; %%vector for allocating # of orange ecDNA copies in each simulation

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
        if xi<prob
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
        if xi<prob
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
        if xi<prob   
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
        if xi<prob
            appob4=appob4+1;
            appor4=appor4-1;
        end
    end
    end
    Ryellow=[Ryellow,appob1+appob3];
    Rred=[Rred,appor1+appor3];
    Ryellow=[Ryellow,appob2+appob4];
    Rred=[Rred,appor2+appor4];
    
    
    Ryellow(j)=[]; %mother removed
    Rred(j)=[];

    else  %ecDNA negative divides
        Ryellow=[Ryellow,0];
        Rred=[Rred,0];
    end

end

%just counting mixed cells
Ryellow2=[];
Rred2=[];

for i=1:length(Ryellow)
    if Ryellow(i)~=0 && Rred(i)~=0
        Ryellow2=[Ryellow2,Ryellow(i)];
        Rred2=[Rred2,Rred(i)];
    end
end

%arranging mixed cells: counting how many cells I have for each configuration
mix=[Rred2' Ryellow2'];
[C,ia,ic]=unique(mix,"rows");
a_counts = accumarray(ic,1);
value_counts = [C, a_counts];
[l,h]=size(value_counts);
for i=1:l
    summa=value_counts(i,1)+value_counts(i,2);
    K=[K,summa];
    %K1=[K1,summa];
    %K=[K,(value_counts(i,1)+value_counts(i,2))/2];
    p=[value_counts(i,1)/summa,value_counts(i,2)/summa]; 
    f=-sum(p.*log(p)); %shannon index formula
    H=[H,f];
end



end

%plot(K,H,'o');

s0=0.1;
Ad=1;
%Plot the density color plot
coldotplot(K,H,s0,Ad);


hold on
l=max(K);
x=linspace(0,l,l);
y=0.6925*ones(1,l);
plot(x,y);  %this is the line of the most frequent shannon index value

appo=[];
appok=[];
for i=2:l
    for j=1:floor(i/2)
        t=[j/i,(i-j)/i];
        f=-sum(t.*log(t));
        appo=[appo,f];
        appok=[appok,i];
    end
end
hold on
h = scatter(appok,appo);  %scatter plot of all the possible shannon indices
axis([0 40.5 0.09 0.71]);
xlabel("Number of ecDNA copies");
ylabel("Shannon index for every cell");



function coldotplot(x,y,s0,Ad)
% Color Scatter Plot for random data point visualization. It mimics a
% continuous 2D probability distribution.
% coldotplot(x,y,s0,Ad) creates at scatterplot with dots of sizes that
% correspond to their density in the swarm of points. The larger dots will
% also have a more "hot" color in the dense particle region.
% The data x and y are vectors of the same size, s0 is a parameter of local
% radii around each datapoint (defalult = 0.5). Ad is a visualization
% parameter for the area of the weighted dots (default = 1).
% Warning: May be slow for very large sizes of x and y.
%
% % Example:
% N=1000;
% x=randn(1,N);
% y=10*randn(1,N);
% s0=0.5;
% Ad=0.2;
% % Plot the density color plot
% coldotplot(x,y,s0,Ad)
%
% Created by: Per Sundqvist, ABB Corporate Research (Swe), 2010-12-07

radiix=s0*std(x);
radiiy=s0*std(y);

phi=linspace(0,2*pi,8);
xv0=radiix*cos(phi);
yv0=radiiy*sin(phi);

colval=zeros(size(x));
for j=1:length(x)
    xv=xv0+x(j);
    yv=yv0+y(j);
    in = inpolygon(x,y,xv,yv);
    colval(j)=sum(in);
end

[dummy,ixcol]=sort(colval);
figure, scatter(x(ixcol),y(ixcol),40,colval(ixcol),'filled')
colormap
colormap copper;
axis tight;
%set(gca,'Color',[0.6784 0.9216 1],'FontSize',20)
end