clf
clear
filename='SEratio1000secs';
SAVE=0;
fontsize=15;
fontname='arial';
markersize=8;
markeredge=markersize/8;

vsub=[14];%subject vein

route=6;%target parameter
% route==1: Area         
% route==2: Perimeter
% route==3: Total LD
% route==4: AP ratio
% route==5: RAP ratio
% route==6: A/sigLd

tageM=7;%target membrane
normalize=2;%1:with actual, 2:maximum
percentage=0;%display as percentage?
compare=1;%compareM6 and 7. 0:no 1:multiply 2:ratio(6/7) 3:multiply and sqrt
rnd=0;%check if the final plot based on the rounded values
rndgt=5;%round in
decimal=1;%if decimal value you are going to show
climmax=1.0;%maximum value of the color bar
climmin=0.0;%minimum value of the color bar

Numcontrs=20;%number of contours drawn in the colormap

cticknum=5;
tick=(0:0.1:1);%how frequently tick appears
sections=1000;

orion=1;
maxon=1;
minon=0;
labelon=1;
climon=1;

%Define parameters---------------------------------------------------------
Parameters=importdata('VeinDimensions.xlsx');

Membrane=importdata('Membrane_Complete_Fly.xlsx');
Junc=importdata('Junctions_Complete_Fly.xlsx');
Direction=importdata('Direction_Complete_Fly.xlsx');
VCircuits=importdata('VCircuits_Complete_Fly.xlsx');

Numv=Parameters.data(:,1);%Vein Numbers
Do=Parameters.data(:,2);%(um)Inner diameters of the wing veins
Lo=Parameters.data(:,3);%(um)Lengths of the wing veins
Ldo=Lo.*Do;

Numm=Membrane.data(:,1);
Amo=Membrane.data(:,2);

Vnum=numel(Numv);
Junctions=Junc.Vcom;
JunctionsS1=Junc.Direction;

%preparation for V14 shift-------------------------------------------------

LoT=Lo;%Keep original value in Lo
Ld=Ldo;%Keep original value in Lo
AmoT=Amo;%The original membrane areas in the triangle

%constants
LaoT=LoT(22)+LoT(15)+LoT(12);
LboT=LoT(13)+LoT(16);

degree=18;
radian=deg2rad(degree);
sine=sin(radian);
cosine=cos(radian);

LoT(6)=sqrt((sine*LaoT)^2+(cosine*LaoT-LboT)^2);
LcoT=LoT(6)+LaoT+LboT;

AmoT67=0.5*LaoT*LboT*sine;

so=LoT(15)/(LaoT-LoT(22));
to=LoT(16)/LboT;
xo=LoT(15)+LoT(22);
yo=LoT(16);
zo=sqrt(xo^2+yo^2-2*xo*yo*cosine);
LoT(14)=zo;

AmoT(7)=0.5*xo*yo*sine;
AmoT(6)=AmoT67-AmoT(7);

LdoT=LoT.*Do;%Keep original value in Lo

imax=numel(VCircuits(:,1));%The number of rows
jmax=numel(VCircuits(1,:));%The number of columns
PerioT=zeros(imax,1);%A box for total kQi in each circuit
Perio=zeros(imax,1);%A box for total kQi in each circuit

AsuploT=zeros(imax,1);%A box for total kQi in each circuit
Asuplo=zeros(imax,1);%A box for total kQi in each circuit
%% 

for i=1:imax
    for j=1:jmax
        if ~isnan(VCircuits(i,j))
            PerioT(i,1)=PerioT(i,1)+LoT(VCircuits(i,j),1);
            Perio(i,1)=Perio(i,1)+Lo(VCircuits(i,j),1);
            AsuploT(i,1)=AsuploT(i,1)+LdoT(VCircuits(i,j),1);
            Asuplo(i,1)=Asuplo(i,1)+Ldo(VCircuits(i,j),1);
        end
    end 
end
RAPo=sqrt(Amo)./Perio;
RAPoT=sqrt(AmoT)./PerioT;
APo=Amo./Perio;
APoT=AmoT./PerioT;
Aevao=2.*Amo;
AevaoT=2.*AmoT;
Rseo=Aevao./Asuplo;
RseoT=AevaoT./AsuploT;

transAm=zeros(sections-2,sections-2,numel(Amo));
transPeri=zeros(sections-2,sections-2,numel(Amo));
transsigLd=zeros(sections-2,sections-2,numel(Amo));

sbox=(1:sections-1)/sections;
tbox=(1:sections-1)/sections;

L=LoT;
Am=AmoT;
%preparation for V14 shift-------------------------------------------------



for s=1:sections-1
    for t=1:sections-1
        L(15)=(LaoT-LoT(22))*s/sections;
        L(16)=LboT*t/sections;
        
        x=L(15)+L(22);
        y=L(16);
        z=sqrt(x^2+y^2-2*x*y*cosine);
        
        L(12)=LaoT-x;
        L(13)=LboT-y;
        L(14)=z;
        Ld=L.*Do;
        
        Am(7)=0.5*x*y*sine;
        Am(6)=AmoT67-Am(7);
        
        Peri=zeros(imax,1);%A box for total kQi in each circuit
        sigLd=zeros(imax,1);%A box for total kQi in each circuit

        for i=1:imax
            for j=1:jmax
                if ~isnan(VCircuits(i,j))
                    Peri(i,1)=Peri(i,1)+L(VCircuits(i,j),1);
                    sigLd(i,1)=sigLd(i,1)+Ld(VCircuits(i,j),1);
                end
            end 
        end
        
        transAm(s,t,:)=Am(:);
        transPeri(s,t,:)=Peri(:);
        transsigLd(s,t,:)=sigLd(:);
    end
end

transRAP=sqrt(transAm)./transPeri;
transAP=transAm./transPeri;
transAsigLd=2*transAm./transsigLd;


   if route==1
        tageparaBox=transAm;
        tageoT=AmoT; 
        
   elseif route==2
        tageparaBox=transPeri;
        tageoT=PerioT;
   
   elseif route==3
        tageparaBox=transsigLd;
        tageoT=AsuploT;
   
   elseif route==4
        tageparaBox=transAP;
        tageoT=APoT;
   
   elseif route==5
        tageparaBox=transRAP;
        tageoT=RAPoT;
        
   elseif route==6
        tageparaBox=transAsigLd;
        tageoT=RseoT;
         
   end 

target=DefTarget(tageparaBox,tageoT, tageM, normalize, percentage, compare);

[tagemax,premaxindex]=max(target(:));
[tagemin,preminindex]=min(target(:));

[premaxs, premaxt] = ind2sub(size(target), premaxindex);
[premins, premint] = ind2sub(size(target), preminindex);

maxs=premaxs/sections;
maxt=premaxt/sections;
mins=premins/sections;
mint=premint/sections;

rso=round(so*sections);
rto=round(to*sections);

tageori=target(rso, rto);

tagefloor=floor(tagemin);
tageceil=ceil(tagemax);
if rnd==1
    target=round(target,rndgt,"significant");
end

if climon==1
    tageindex=climmin:(climmax-climmin)/Numcontrs:climmax;
    
    if decimal==1
        [C, h]=contourf(tbox, sbox,target,tageindex,"LabelFormat","%0.2f");
    else
        [C, h]=contourf(tbox, sbox,target,tageindex);
    
    end
    climit=[climmin climmax];
    clim(climit)
    
    
 color=copper;
        colorR=color(:,1);
        colorG=color(:,2);
        colorB=color(:,3);
    rmax=1;
    rmin=0.99;
    color(:,1)=(rmin:(rmax-rmin)/255:rmax);
    
    gmax=0.99;
    gmin=0.50;
    color(:,2)=(gmin:(gmax-gmin)/255:gmax);
    
    bmax=0.98;
    bmin=0.15;
    color(:,3)=(bmin:(bmax-bmin)/255:bmax);
    
    color=flipud(color);
    colormap(color)
   
    c=colorbar;
    c.Box='on';
    c.LineWidth=0.5;
    
    ctick=(climmax-climmin)/cticknum;%how frequently tick appears

    c.Ticks=(climmin:ctick:climmax);

    % c.Label.String ='{\it WSI}/{\it WSI}{_o} [%]';
    % c.Label.String ='{\it Q}/{\it Q}{_o} [%]';
    c.FontName=fontname;
    c.FontSize=0.8*fontsize;

    
    %mesh(tbox, sbox, transQn(:,:,NumVein))
    %mesh(tbox, sbox, tQLperAm67n)

    
elseif climon==0
    [C, h]=contourf(tbox, sbox,target);
    color=bone;
    color=flipud(color);
    colormap(color)
    colorbar
    c=colorbar;
    c.Box='on';
    c.LineWidth=0.8;

    % c.Label.String ='{\it WSI}/{\it WSI}{_o} [%]';
    % c.Label.String ='{\it Q}/{\it Q}{_o} [%]';
    c.FontName=fontname;
    c.FontSize=0.8*fontsize;

    
    %mesh(tbox, sbox, transQn(:,:,NumVein))
    %mesh(tbox, sbox, tQLperAm67n)

end

if labelon
    clabel(C,h,'LabelSpacing',200,'fontsize',fontsize*0.6,'fontname',fontname)
end 


if orion==1
    hold on
    plot3(to,so,0,'ro','MarkerSize',markersize,'LineWidth',markeredge,'MarkerEdgeColor',[.15 .0 .0],...
        'MarkerFaceColor',[.65 0.1 .1])
end 

if maxon==1
    
    hold on
    plot3(maxt,maxs,0,'ro','MarkerSize',round(markersize/1.5),'LineWidth',markeredge/1.5,'MarkerEdgeColor',[.05 .05 .05],...
        'MarkerFaceColor',[.3 0.3 .3])
end

if minon==1
    hold on
    plot3(mint,mins,0,'ro','MarkerSize',round(markersize/1.5),'LineWidth',markeredge/1.5,'MarkerEdgeColor',[.05 .05 .05],...
        'MarkerFaceColor',[.9 0.9 .9])
end


xlim([0 1])
ylim([0 1])
xticks(tick)
yticks(tick)
xticklabels({})
yticklabels({})
% zlim([50 150])
set(gca,'FontSize',fontsize*0.8,'FontName',fontname);
% xlabel('{\it t} [%]','FontSize',fontsize,'FontName',fontname)
% ylabel('{\it s} [%]','FontSize',fontsize,'FontName',fontname)
% zlabel('Total pressure difference [Pa]','FontSize',fontsize,'FontName',fontname)
% zlabel('{\it WSI}/{\it WSI}{_o} [%]','FontSize',fontsize,'FontName',fontname)
%zlabel('{\it Q}/{\it Q}{_o} [%]','FontSize',fontsize,'FontName',fontname)

ax=gca;
ax.Position(1)=0.1;
ax.PlotBoxAspectRatio(2)=1;
ax.PlotBoxAspectRatio(3)=1;

if SAVE==1
    save(filename)
end

function target=DefTarget(tageparaBox,tageoT, tageM, normalize, percentage, compare) 
        if normalize==1
            for c=1:numel(tageoT)
                tageparaBox(:,:,c)=tageparaBox(:,:,c)./tageoT(c)*(1+99*percentage);       
            end                      
        end
                     
        if compare==1
            target=tageparaBox(:,:,6).*tageparaBox(:,:,7);
        elseif compare==2
            target=tageparaBox(:,:,6)./tageparaBox(:,:,7);
        elseif compare==3
            target=sqrt(tageparaBox(:,:,6).*tageparaBox(:,:,7));            
        else
            target=tageparaBox(:,:,tageM);
        end
        
        if normalize==2
            target=target./max(target,[],"all");
        end
end
