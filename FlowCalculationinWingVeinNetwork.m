
clf
clear
tic

%--
%Decide do or undo-------------------------------------------
wholerepeat=1;
savevalues=0;
%Decide do or undo-------------------------------------------
%--

%--
%Define calculation conditions---------------------------------------------------------
sections=1000;
NumPlotson1side=sections-1;

TageCVList=["excv","pcv"];
%["excv","pcv"]

%Define calculation conditions---------------------------------------------------------
%--

%--
%Define network parameters-------------------------------------------
species='D.melanogaster';
%D.cyrtoloma
%D.melanogaster


veintype="Normal";
% veintype='NoExCV';
%NoPCVExCV
%NoPCV
%Normal
%NoExCV
%Normal4Matrix
%NoCV_Model
%NoExCV_Model
%ShiftCV4Matrix
%LargeD

%Define network parameters-------------------------------------------
%--
    
FluidicParameters.InflowRate=520;
%abs inflow rate
FluidicParameters.Density=1002*1e-18;
%(kg/um^3)Density of hemolymph;%(kg/um^3)Density of hemolymph
FluidicParameters.Viscosity=0.0013*1e-6;
%(kg/(um*s))viscousity of hemolymph(Pa*s=kg*m/s^2/m^2*s=kg/(m*s)=10^-6 kg/(um*s))
FluidicParameters.GravitationalAcc=9.8*10^6;
%(um/s^2)gravitational acceleration

%%
%-----------------------------------------------
Path4Folder='..\CircuitAnalogy\Hardy-Cross\HCprarameters\Fruitfly\';
path=append(Path4Folder,species,'\',species,'_',veintype,'\');
Atrbt=append('_',species,'_',veintype,'.xlsx');
%create path to treat files or save files

VeinGeometryOriginal=readtable(append(path,'VeinGeo',Atrbt),...
    'VariableNamingRule','preserve');
[VeinNames,~]=SeparateName(VeinGeometryOriginal);
VeinGeometryOriginal.CrossSectionalArea=...
    pi*VeinGeometryOriginal.InnerDiameter.^2/4;
%(um^2)Cross-sectional area of the wing veins

VeinGeometryOriginal.ContactArea=...
    0.5*pi*VeinGeometryOriginal.InnerDiameter.*VeinGeometryOriginal.Length;     
%read and created the table containing the original vein geometry

LoopsOriginal=readtable(append(path,'Loops',Atrbt),...
    'VariableNamingRule','preserve');
LoopsOriginal=mergevars(LoopsOriginal,...
    2:numel(LoopsOriginal(1,:)),'NewVariableName',"Veins");


MemAreaPeriListOriginal=readtable(append(path,'LoopGeoOriginal',Atrbt),...
    'VariableNamingRule','preserve');
[MembraneGeometryOriginal]=ObtainMemGeoOriginal...
    (VeinGeometryOriginal,LoopsOriginal, MemAreaPeriListOriginal);

NodeswithJoinedVeinsOriginal=readtable(append(path,'Nodes',Atrbt));

%grasp data
FlowDirectionatNodesOriginal=readtable(append(path,'NodeDirections',Atrbt),...
    'VariableNamingRule','preserve');

FlowDirectioninLoopsOriginal=readtable(append(path,'LoopDirections',Atrbt),...
    'VariableNamingRule','preserve');
FlowDirectioninLoopsOriginal=mergevars(FlowDirectioninLoopsOriginal,...
    2:numel(FlowDirectioninLoopsOriginal(1,:)),'NewVariableName',"Direction");

[FlowCalcResultOriginal]=FlowRateCalculation...
    (VeinGeometryOriginal,LoopsOriginal,FlowDirectioninLoopsOriginal,...
    FluidicParameters);

%%
if contains(veintype,"Shift")==1
    NodesinVeinFramework=readtable(append...
        (path,'NodesinVeinFramework',Atrbt),'NumHeaderLines',0);
    
    NodesinVeinFramework=mergevars(NodesinVeinFramework,...
        2:numel(NodesinVeinFramework(1,:)),'NewVariableName',"Nodes");
    %read the table shows relationship between general names of veinframes and nodes
    %connecting these frames
    
    NodePositionListOriginal=readtable(append...
        (path,'NodeCoordinates',Atrbt),'NumHeaderLines',0);

    NodePositionListOriginal=mergevars(NodePositionListOriginal,...
        2:numel(NodePositionListOriginal(1,:)),'NewVariableName',"Coordinates");
    %read the table includes the positions of landmark nodes in the simplified model
    
    %%
    [NodePositionListInitial]=InitializeNodePositions...
        (NodePositionListOriginal,TageCVList,VeinNames, LoopsOriginal);
    %%Initialize all node positions including CV nodes
    
    [VeinFrameworkGeometry]=ObtainVeinFrameworkGeometry...
        (NodePositionListInitial,NodesinVeinFramework,...
        TageCVList, LoopsOriginal);
    %%Obtain Geometry of vein framework 
    
    [LoopsInitial, FlowDirectioninLoopsInitial]=ObtainLoopVeinsandFlowDir...
        (NodePositionListInitial,LoopsOriginal, FlowDirectioninLoopsOriginal);
    
    if any(contains(VeinNames,"m"))
        [NodeswithJoinedVeinsInitial,FlowDirectionatNodesInitial]=ObtainJoinedVeinsandFlowDir...
            (NodePositionListInitial,NodeswithJoinedVeinsOriginal,FlowDirectionatNodesOriginal);
    else
        NodeswithJoinedVeinsInitial=NodeswithJoinedVeinsOriginal;
    end
    %Obtain Veins joined at each node reflecting node positions of CVs
    

    [VeinListwithNodesInitial]=ObtainNodesofVein(NodeswithJoinedVeinsInitial,VeinNames);
    %Obtain Nodes forming respective veins from the Joined Vein List
    
    [VeinListwithNodesExtracted]=ExtractNodesofVeinList(VeinListwithNodesInitial);
    %Extract the node pair lists for the veins that have variable lengths
    
    [VeinGeometryInitial,Origins]=ObtainVeinGeometry...
        (VeinFrameworkGeometry,VeinGeometryOriginal,NodePositionListInitial,...
        VeinListwithNodesExtracted,TageCVList,LoopsInitial);
    %obtain the vein lengths defined by the initial node positions
    
    
    MemAreaPeriListInitial=readtable(append(path,'LoopGeoInitial',Atrbt));
    [MembraneGeometryPreInitial]=ObtainMemGeoOriginal...
        (VeinGeometryInitial,LoopsInitial, MemAreaPeriListInitial);
    
    [MembraneGeometryInitial]=ObtainMembraneGeometry...
        (MembraneGeometryPreInitial,VeinGeometryInitial,LoopsInitial,...
        NodePositionListInitial, NodesinVeinFramework,TageCVList);
        
    [FlowCalcResultInitial]=FlowRateCalculation...
        (VeinGeometryInitial,LoopsInitial,FlowDirectioninLoopsInitial,...
        FluidicParameters);
    FlowCalcResultInitial.Qi_Qo=...
        FlowCalcResultInitial.FlowRate./FlowCalcResultOriginal.FlowRate;
    FlowCalcResultInitial.pi_po=...
        FlowCalcResultInitial.PressureLoss_Pa./FlowCalcResultOriginal.PressureLoss_Pa;
    
%     plot(NodePositionListInitial.Coordinates(:,1),NodePositionListInitial.Coordinates(:,2),...
%         MarkerFaceColor="k",MarkerEdgeColor="k",Marker="o",MarkerSize=5)
%     hold on
%     
%     xlim([0 2735])
%     ylim([0 2325])
%     pbaspect([2735 2325 1])
%     hold on
    
    
    %%
    
    if wholerepeat==1
        NumTageCVs=numel(TageCVList);
        NumVeins=numel(VeinNames);
        NumLoops=numel(LoopsInitial.LoopName);
        for tagecv=1:NumTageCVs
            TageCV=TageCVList(tagecv);
            

            [NumVcp]=FindTargetNumofVe(LoopsInitial,TageCV);
      
            transtotalp=zeros(sections-1,sections-1);
            transl=zeros(sections-1,sections-1,NumVeins);
            transQ=zeros(sections-1,sections-1,NumVeins);
            transp=zeros(sections-1,sections-1,NumVeins);
            transContactArea=zeros(sections-1,sections-1,NumVeins);
    
            transAm=zeros(sections-1,sections-1,NumLoops);
            transAsupply=zeros(sections-1,sections-1,NumLoops);
            transPerimeter=zeros(sections-1,sections-1,NumLoops);
            
           
            parfor sectionP=1:NumPlotson1side
                for sectionA=1:NumPlotson1side
                    [NodePositionListShift]=ObtainCVNodeCoordinates...
                    (NodePositionListInitial,VeinFrameworkGeometry,LoopsInitial,...
                    sectionA,sectionP,sections, TageCV);
                    
                    [LoopsShift, FlowDirectioninLoopsShift]=ObtainLoopVeinsandFlowDir...
                        (NodePositionListShift,LoopsInitial, FlowDirectioninLoopsInitial);
                   
                    [NodeswithJoinedVeinsShift,FlowDirectionatNodesShift]=...
                        ObtainJoinedVeinsandFlowDir...
                        (NodePositionListShift,NodeswithJoinedVeinsInitial,FlowDirectionatNodesInitial);
                
                    
                    %Obtain Veins joined at each node reflecting node positions of CVs
                    
                    [VeinListwithNodesShift]=ObtainNodesofVein(NodeswithJoinedVeinsShift,VeinNames);
                    %Obtain Nodes forming respective veins from the Joined Vein List
                    
                    [VeinListwithNodesExtracted]=ExtractNodesofVeinList(VeinListwithNodesShift);
                    %Extract the node pair lists for the veins that have variable lengths
                    
                    [VeinGeometryShift]=ObtainVeinGeometry...
                        (VeinFrameworkGeometry,VeinGeometryInitial,NodePositionListShift,...
                        VeinListwithNodesExtracted,TageCVList,LoopsShift);
                    %obtain the vein lengths defined by the initial node positions
                    
                    [FlowCalcResultShift]=FlowRateCalculation...
                        (VeinGeometryShift,LoopsShift,FlowDirectioninLoopsShift,...
                        FluidicParameters);
    
                    transQ(sectionA,sectionP,:)=FlowCalcResultShift.FlowRate;
                    transp(sectionA,sectionP,:)=FlowCalcResultShift.PressureLoss_Pa;
                    transtotalp(sectionA,sectionP)=FlowCalcResultShift.TotalPressureLoss_Pa(1);
                    transl(sectionA,sectionP,:)=FlowCalcResultShift.Length;
    
                    l_shift=FlowCalcResultShift.Length;
                    d_shift=FlowCalcResultShift.InnerDiameter;
                    ContactArea_shift=0.5*pi*l_shift.*d_shift;
                    transContactArea(sectionA,sectionP,:)=ContactArea_shift;
                                    
                    [MembraneGeometryShift]=ObtainMembraneGeometry...
                        (MembraneGeometryInitial,VeinGeometryShift,LoopsShift,...
                        NodePositionListShift, NodesinVeinFramework,TageCVList);
              
                    transAm(sectionA,sectionP,:)=MembraneGeometryShift.Area;
                    transPerimeter(sectionA,sectionP,:)=MembraneGeometryShift.Perimeter;
                    transAsupply(sectionA,sectionP,:)=MembraneGeometryShift.SupplyArea;
                    
                    "Target Vein："+TageCV+newline+...
                        "Node A："+num2str(sectionA)+"/"+num2str(NumPlotson1side)+...
                        newline+"Node P："+num2str(sectionP)+"/"+num2str(NumPlotson1side)
                
                end
                
                
            end

                    transAevapo=transAm*2;
                    transRAP=sqrt(transAm)./transPerimeter;
                    transRse=transAevapo./transAsupply;
                    LoopNums=(1:NumLoops);
                    LoopNume=LoopNums(contains(LoopsInitial.LoopName,append(num2str(NumVcp),"e")));
                    LoopNumb=LoopNums(contains(LoopsInitial.LoopName,append(num2str(NumVcp),"b")));
                    transsqrtRse=sqrt(transRse(:,:,LoopNume).*transRse(:,:,LoopNumb));
                    transsqrtRAP=sqrt(transRAP(:,:,LoopNume).*transRAP(:,:,LoopNumb));
                     
                    save (append(species,"_",veintype,"_",TageCV,"_",num2str(sections),"sections.mat"));
                    writetable(FlowCalcResultInitial,append(species,"_FlowCalcResult_Initial",veintype,".xlsx"));
                    writetable(MembraneGeometryInitial,append(species,"_MembraneGeometry_Initial",veintype,".xlsx"));

        end
    
    end

    % poolobj = gcp('nocreate');
    % delete(poolobj);


end

if savevalues==1
    save (append(species,"_",veintype,".mat"))
    writetable(FlowCalcResultOriginal,append(species,"_FlowCalcResult_",veintype,".xlsx"));
    writetable(MembraneGeometryOriginal,append(species,"_MembraneGeometry_",veintype,".xlsx"));
end
toc
    
%% 

%Definitions of functions--------------------------------------------------
%---------------------------------------------------------------------------

%%
function [SpecificNames]=SpecializeNames(GeneralNames,NumVcp)
    SpecificNames=GeneralNames;
    
    SpecificNames=strrep(SpecificNames,"Ve",append("Ve",num2str(NumVcp)));
    SpecificNames=strrep(SpecificNames,"Vb",append("Vb",num2str(NumVcp)));

    SpecificNames=strrep(SpecificNames,"a",num2str(NumVcp-1));
    SpecificNames=strrep(SpecificNames,"p",num2str(NumVcp));
end


%--
%Separates Names from Contents of lists--------------------------------------------
function [Names, Contents]=SeparateName(ParaTable)
    Names=ParaTable(:,1);
    Names=table2array(Names);
    if class(Names)=="cell"
        Names=string(Names);
    end
    
    Contents=ParaTable(:,2:end);
    Contents=table2array(Contents);
    
    if class(Contents)=="cell"
        Contents=string(Contents);
    end
end
%Separates Names from Contents of lists--------------------------------------------
%--
%%
function [VeNum]=FindTargetNumofVe(Loops, TageCV)
    [~,LoopVeinComs]=SeparateName(Loops);
    %translate table to matrix and separate names from matrix
    ContainTageVein=sum(LoopVeinComs==TageCV,2);
    %find membrane cells containing the target vein
    SrndVeins=LoopVeinComs(ContainTageVein==1,:);
    %find the veins surrounding the cell
   
    Ve=SrndVeins(contains(SrndVeins,"Ve"));
    %extract name of edge and base veins in the cell region
    pat = digitsPattern;
    VeNum = str2double(extract(Ve,pat));
end

%%
%--
%Translate vein names to numbers in vein connection lists--------------------------------------------
function [ComponentsasNums]=TranslateName2Num(ComponentsasNames,ComponentNames)
    NumTargets=numel(ComponentsasNames(:,1));%the number of Nodes
    MaxNumofComponents=numel(ComponentsasNames(1,:));%The maximum numbers connected at node
    ComponentsasNums=NaN(NumTargets, MaxNumofComponents);%nan box for translate names to numbers
    for i=1:NumTargets
        for j=1:MaxNumofComponents
            if ComponentsasNames(i,j)~=""
                f=find(ComponentNames==ComponentsasNames(i,j));
                ComponentsasNums(i,j)=f;
            end
        end 
    end
end
%Translate vein names to numbers in vein connection lists--------------------------------------------
%--
%%
function [NodePositionListShift]=ObtainCVNodeCoordinates...
    (NodePositionListInitial,VeinFrameworkGeometry,Loops,...
    sectionA,sectionP,sections, TageCV)
    
    NodePositionListShift=NodePositionListInitial;
    
    FrameNames=VeinFrameworkGeometry.VeinFrame;
    FrameLength=VeinFrameworkGeometry.Length;

    CoordinateList=NodePositionListShift.Coordinates;
    NodeNameList=NodePositionListShift.Node;
    
    [NumVcp]=FindTargetNumofVe(Loops, TageCV);
    NumVca=NumVcp-1;
    sectionAP=([sectionA, sectionP]);
    signs=(["a","p"]);
    
    
    for NumTageVc=NumVca:NumVcp%anterior Vc first, and posterior Vc later
        NumAorP=NumTageVc-NumVcp+2;%1:A,2:P
        AorP=signs(NumAorP);
        section=sectionAP(NumAorP);
        FractionSectioned=section/sections;

        NameTageFrameEdge=append("Vc",num2str(NumTageVc),"E");
        NameTageFrameBase=append("Vc",num2str(NumTageVc),"B");
        NcName=append("Nc",num2str(NumTageVc));
        NbName=append("Nb",num2str(NumTageVc));
        NeName=append("Ne",num2str(NumTageVc));
        NcvName=append("N",TageCV,AorP);

        lcE=FrameLength(FrameNames==NameTageFrameEdge);
        lcB=FrameLength(FrameNames==NameTageFrameBase);
        Lc=lcE+lcB;

        lcb=Lc*FractionSectioned;
        
        Nc=CoordinateList(NodeNameList==NcName,:);
        Nb=CoordinateList(NodeNameList==NbName,:);
        Ne=CoordinateList(NodeNameList==NeName,:);
      
        if lcb<=lcB
            ratio=lcb/lcB;
            Ncv=(Nc-Nb)*ratio+Nb;
        else
            ratio=(lcb-lcB)/lcE;
            Ncv=(Ne-Nc)*ratio+Nc;
        end
        NodePositionListShift.Coordinates((NodeNameList==NcvName),:)=Ncv;
    end 
end


%%
function [VeinGeometry,Origins]=ObtainVeinGeometry...
    (VeinFrameworkGeometry,VeinGeometryOriginal,NodePositionList,...
    VeinListwithNodesExtracted,...
    TageCVList,Loops)

    VeinGeometry=VeinGeometryOriginal;
    VeinNames=string(VeinGeometry.Vein);
    Origin=zeros(numel(TageCVList),1);
    
    for tagecv=1:numel(TageCVList)
        TageCV=TageCVList(tagecv);
        [NumVcp]=FindTargetNumofVe(Loops, TageCV);
        
        CoordinateList=NodePositionList.Coordinates;
        NodeNameList=NodePositionList.Node;
    
        for i=1:numel(VeinListwithNodesExtracted.Vein)
            nodepair=VeinListwithNodesExtracted.Nodes(i,:);
            Node1=nodepair(1);
            Node2=nodepair(2);
            Coordinates1=CoordinateList(NodeNameList==Node1,:);
            Coordinates2=CoordinateList(NodeNameList==Node2,:);
            VeinListwithNodesExtracted.Length(i)=norm(Coordinates2-Coordinates1);
        end
        %calculate lengths of veins in frame as straight lines connecting corresponding nodes 
        
        VeinListwithNodesExtractedRedefined=VeinListwithNodesExtracted;
        
        Vein=VeinListwithNodesExtractedRedefined.Vein;
        Nodes=VeinListwithNodesExtractedRedefined.Nodes;
          
        FrameNames=VeinFrameworkGeometry.VeinFrame;
        FrameLength=VeinFrameworkGeometry.Length;
        FramesIntersect=(["";""]);
    
        
        for NumTageVc=NumVcp-1:NumVcp%anterior Vc first, and posterior Vc later
            
            NameTageFrameEdge=append("Vc",num2str(NumTageVc),"E");
            NameTageFrameBase=append("Vc",num2str(NumTageVc),"B");
            lcE=FrameLength(FrameNames==NameTageFrameEdge);
            lcB=FrameLength(FrameNames==NameTageFrameBase);
            Lc=lcE+lcB;
            

            NameTageVeinEdge=append("Vc",num2str(NumTageVc),"e");
            NameTageVeinBase=append("Vc",num2str(NumTageVc),"b");
            
            nodepaire=Nodes((Vein==NameTageVeinEdge),:);
            nodepairb=Nodes((Vein==NameTageVeinBase),:);
        
            NcveName=nodepaire(contains(nodepaire,"cv"));
            NcvbName=nodepairb(contains(nodepairb,"cv"));
        
            NcName=append("Nc",num2str(NumTageVc));
            Nc=CoordinateList(NodeNameList==NcName,:);
            Ncve=CoordinateList(NodeNameList==NcveName,:);
            Ncvb=CoordinateList(NodeNameList==NcvbName,:);
        
            lce=lcE+norm(Ncve-Nc)*sign(Nc(1)-Ncve(1));
            lcb=lcB+norm(Ncvb-Nc)*sign(Ncvb(1)-Nc(1));
            
            VeinListwithNodesExtractedRedefined.Length...
                (VeinListwithNodesExtractedRedefined.Vein==NameTageVeinEdge)=lce;
            VeinListwithNodesExtractedRedefined.Length...
                (VeinListwithNodesExtractedRedefined.Vein==NameTageVeinBase)=lcb;
        
            if NumTageVc==5
                lcm=Lc-(lce+lcb);
                VeinListwithNodesExtractedRedefined.Length...
                (VeinListwithNodesExtractedRedefined.Vein=="Vc5m")=lcm;

            end
            

            if lcb<=lcB
                FramesIntersect(NumTageVc-NumVcp+2,1)=NameTageFrameBase;
            else 
                FramesIntersect(NumTageVc-NumVcp+2,1)=NameTageFrameEdge;
            end

            Origin(tagecv,NumTageVc-NumVcp+2)=lcB/Lc;
        end
        
        
        Ncva=CoordinateList(NodeNameList==append("N",TageCV,"a"),:);
        Ncvp=CoordinateList(NodeNameList==append("N",TageCV,"p"),:);
        Nint=ObtainIntersection(FramesIntersect,VeinFrameworkGeometry, NodePositionList);
        if (Ncva(1)<Nint(1))||(Ncvp(1)<Nint(1))
            lcv=NaN;
        else
            lcv=norm(Ncva-Ncvp);
        end
        
        VeinListwithNodesExtractedRedefined.Length...
                (VeinListwithNodesExtractedRedefined.Vein==TageCV)=lcv;
        Velist=VeinFrameworkGeometry(contains(VeinFrameworkGeometry.VeinFrame,"Ve"),:);
        Vblist=VeinFrameworkGeometry(contains(VeinFrameworkGeometry.VeinFrame,"Vb"),:);
        VeVblist=([Velist;Vblist]);
        VeVblist=renamevars(VeVblist,"VeinFrame","Vein");
        VeinListwithNodesExtractedRedefined=([VeinListwithNodesExtractedRedefined;VeVblist]);
    
        vlengthExtracted=VeinListwithNodesExtractedRedefined.Length;
        vnamesExtracted=VeinListwithNodesExtractedRedefined.Vein;
        
        for scan=1:numel(VeinListwithNodesExtractedRedefined.Vein)
            VeinGeometry.Length(VeinNames==vnamesExtracted(scan))...
                =vlengthExtracted(scan);
        end
        
        %Judge Vc5m diameter
        di=VeinGeometry.InnerDiameter;
        % VeinGeometryInitial.InnerDiameter(VeinNames=="Vc5m")=mean([di(VeinNames=="Vc5e"),di(VeinNames=="Vc5b")]);
        %constant diameter
        if any(contains(VeinNames,"m"))
            x_excvp=CoordinateList(NodeNameList=="Nexcvp",1);
            x_pcva=CoordinateList(NodeNameList=="Npcva",1);
            x_c=CoordinateList(NodeNameList=="Nc5",1);
            x_basal=min([x_pcva,x_excvp]);
            if x_basal<x_c
                VeinGeometry.InnerDiameter(VeinNames=="Vc5m")=di(VeinNames=="Vc5b");
            else
                VeinGeometry.InnerDiameter(VeinNames=="Vc5m")=di(VeinNames=="Vc5e");
            end
            %Variable diameter based on the position of the most basal CV node
        end

        VeinGeometry.ContactArea=0.5*pi*...
        VeinGeometry.InnerDiameter.*VeinGeometry.Length;     

%         plot(Ncve(1),Ncve(2),...
%         MarkerFaceColor="b",MarkerEdgeColor="k",Marker="o",MarkerSize=5)
%         hold on
%         
%         plot(Ncvb(1),Ncvb(2),...
%         MarkerFaceColor="r",MarkerEdgeColor="k",Marker="o",MarkerSize=5)
%         
%         xlim([0 2735])
%         ylim([0 2325])
%         pbaspect([2735 2325 1])
%         hold on    

    end
    CrossVeins=TageCVList.';
    Origins=table(CrossVeins,Origin);
end
%%
function [NodeswithJoinedVeinsNew, FlowDirectionsatNodesNew]=ObtainJoinedVeinsandFlowDir...
    (NodePositionList,NodeswithJoinedVeinsPrev, FlowDirectionsatNodesPrev)
    

        CoordinateList=NodePositionList.Coordinates;
        x_excvp=CoordinateList(NodePositionList.Node=="Nexcvp",1);
        x_pcva=CoordinateList(NodePositionList.Node=="Npcva",1);
    
        NodeswithJoinedVeinsNew=NodeswithJoinedVeinsPrev;
        FlowDirectionsatNodesNew=FlowDirectionsatNodesPrev;
        
        [NodeNames,VeinComponents]=SeparateName(NodeswithJoinedVeinsNew);
        contains_m=contains(VeinComponents,"m");
        
        if any(contains_m,"all")
            if x_pcva<=x_excvp
                vein4excvp=(["excv","Vc5e", "Vc5m"]);
                vein4pcva=(["Vc5m","pcv","Vc5b"]);
        
                Dir4excvp=([1,1,-1]);
                Dir4pcva=([1,-1,-1]);
            else
                vein4excvp=(["excv" "Vc5m" "Vc5b"]);
                vein4pcva=(["Vc5e","pcv","Vc5m"]);
        
                Dir4excvp=([1,1,-1]);
                Dir4pcva=([1,-1,-1]);
            end
            
            for i=1:numel(vein4pcva)
                NodeswithJoinedVeinsNew(NodeNames=="Nexcvp",i+1)=cellstr(vein4excvp(1,i));
                NodeswithJoinedVeinsNew(NodeNames=="Npcva",i+1)=cellstr(vein4pcva(1,i));
        
                FlowDirectionsatNodesNew(NodeNames=="Nexcvp",i+1)=table(Dir4excvp(1,i));
                FlowDirectionsatNodesNew(NodeNames=="Npcva",i+1)=table(Dir4pcva(1,i));
            end

        end
  
end

%%
function [LoopsNew, FlowDirectionsinLoopsNew]=ObtainLoopVeinsandFlowDir...
    (NodePositionList,LoopsPrev, FlowDirectionsinLoopsPrev)

    LoopsNew=LoopsPrev;
    FlowDirectionsinLoopsNew=FlowDirectionsinLoopsPrev;
    [LoopNames,VeinComponents]=SeparateName(LoopsNew);
    Direction=FlowDirectionsinLoopsNew.Direction;
    contains_m=contains(VeinComponents,"m");

    if any(contains_m,"all")

        CoordinateList=NodePositionList.Coordinates;
        x_excvp=CoordinateList(NodePositionList.Node=="Nexcvp",1);
        x_pcva=CoordinateList(NodePositionList.Node=="Npcva",1);
     
        if x_pcva<=x_excvp
            ContainTageLoopVc5m=contains(LoopNames,(["L5b","L6e"]));
            ContainTageLoopOther=contains(LoopNames,(["L5e","L6b"]));
        else
            ContainTageLoopVc5m=contains(LoopNames,(["L5e","L6b"]));
            ContainTageLoopOther=contains(LoopNames,(["L5b","L6e"]));
    
        end
        
        VeinComponents(ContainTageLoopVc5m,5)=cellstr(["Vc5m" "Vc5m"]);
        VeinComponents(ContainTageLoopOther,5)=cellstr(["" ""]);
        
        Direction(ContainTageLoopVc5m,5)=([1,-1]);
        Direction(ContainTageLoopOther,5)=([nan,nan]);
        
        LoopsNew.Veins=VeinComponents;
        FlowDirectionsinLoopsNew.Direction=(Direction);
    end
end
%%
function [NodesofVeinList]=ObtainNodesofVein(Nodes,VeinNames)
    [NodeNames,VeinsJoined]=SeparateName(Nodes);
    Vein=strings(numel(VeinNames),1);
    Nodes=strings(numel(VeinNames),2);
    
    for scan=1:numel(VeinNames)
        vname=VeinNames(scan);
        RowsContaintheVein=logical(sum(VeinsJoined==vname,2));
        nodepair=NodeNames(RowsContaintheVein);
        Vein(scan)=vname;
        Nodes(scan,:)=nodepair(:);
    end
    NodesofVeinList=table(Vein,Nodes);
end

%%
function [NodesofVariableVeinsExtracted]=ExtractNodesofVeinList(NodesofVeinList)
    Vein=NodesofVeinList.Vein;
    RowsContainVeinName4=strlength(Vein)==4;
    RowsContainPCV=Vein=="pcv";
    RowsContainTageVein=logical(RowsContainPCV+RowsContainVeinName4);
    NodesofVariableVeinsExtracted=NodesofVeinList(RowsContainTageVein,:);
end


%%
function [NodePositionListInitial]=InitializeNodePositions(NodePositionListOriginal,TageCVList,VeinNames, Loops)
   
    Node=strings(2*numel(TageCVList),1);
    Coordinates=zeros(2*numel(TageCVList),2);
    
    for i=1:numel(TageCVList)
        TageCV=TageCVList(i);
        %extract target vein
        
        [NumMemTarget]=FindTargetNumofVe(Loops, TageCV);
        
        NumVcp=NumMemTarget;
        NumVca=NumVcp-1;

        Node(2*i-1,1)=append("N",TageCV,"a");
        Node(2*i,1)=append("N",TageCV,"p");
        Coordinates(2*i-1,:)=NodePositionListOriginal.Coordinates(contains(NodePositionListOriginal.Node,append("c",num2str(NumVca))),:);
        Coordinates(2*i,:)=NodePositionListOriginal.Coordinates(contains(NodePositionListOriginal.Node,append("c",num2str(NumVcp))),:);
    end
    NodePositionListInitialCV=array2table(Node);
    NodePositionListInitialCV.Coordinates=Coordinates;
    NodePositionListInitial=cat(1,NodePositionListOriginal,NodePositionListInitialCV);

end


%%
function [ExtractedNodePositionList]=ExtractTargetNodes(OriginalNodePositionList, TageCV, NumVcp)
    [OriginalNodeNames,~]=SeparateName(OriginalNodePositionList);
    withb=contains(OriginalNodeNames,"b");
    withc=contains(OriginalNodeNames,"c");
    withe=contains(OriginalNodeNames,"e");
    withA=contains(OriginalNodeNames,num2str(NumVcp-1));
    withP=contains(OriginalNodeNames,num2str(NumVcp));
    withCV=contains(OriginalNodeNames,TageCV);
        
    TargetNodeRows=logical(withA.*(withb+withc+withe)+withP+withCV);
    ExtractedNodePositionList=OriginalNodePositionList(TargetNodeRows,["Node","Coordinates"]);
end

function [NodePositionListGeneralized]=GeneralizeNodeNames(ExtractedNodePositionList, TageCV, NumVcp)
    NodePositionListGeneralized=ExtractedNodePositionList;
    [NodeNames,~]=SeparateName(NodePositionListGeneralized);
    NodeNames = strrep(NodeNames,num2str(NumVcp-1),"a");
    NodeNames = strrep(NodeNames,num2str(NumVcp),"p");
    NodeNames(contains(NodeNames,"C"))="C";
    NodeNames(contains(NodeNames,append(TageCV,"a")))="Na";
    NodeNames(contains(NodeNames,append(TageCV,"p")))="Np";
    NodePositionListGeneralized.Node=NodeNames;

end 


%%
function [VeinFramework]=ObtainVeinFramework(NodesinVeinFramework,NodePositionListGeneralInitial)
    [NodeNames, NodeCoordinates]=SeparateName(NodePositionListGeneralInitial);
    [FrameNames, NodesofFrame]=SeparateName(NodesinVeinFramework);
    NumFrames=numel(NodesinVeinFramework(:,1));
    FrameLength=zeros(NumFrames,1);

    for i=1:NumFrames
        N=NodesofFrame(i,:);
        NodePositionBox=NodeCoordinates(contains(NodeNames,N),:);
        FrameLength(i)=norm(NodePositionBox(1,:)-NodePositionBox(2,:));
    end
     
    Ra=FrameLength(contains(FrameNames,"Ra"));
    Rp=FrameLength(contains(FrameNames,"Rp"));
    R=0.5*(Ra+Rp);
    
    Nea=NodeCoordinates(contains(NodeNames,"Nea"),:);
    Nep=NodeCoordinates(contains(NodeNames,"Nep"),:);
    Lchrd=norm(Nea-Nep);
    costheta=(R^2+R^2-Lchrd^2)/(2*R*R);
    theta=acos(costheta);
    %theta_deg=rad2deg(theta);
    le=R*theta;
    VeinFrame=(["Vb";"VcaE";"VcaB";"VcpE";"VcpB"]);
    Length=FrameLength(contains(FrameNames,VeinFrame));
    Nodes=NodesofFrame(contains(FrameNames,VeinFrame),:);
    VeinFrame=(["Ve";VeinFrame]);
    Length=([le;Length]);
    Nodes=([NaN,NaN;Nodes]);
    VeinFramework=table(VeinFrame,Nodes,Length);

end


%%
function [VeinFrameworkGeometry]=ObtainVeinFrameworkGeometry...
    (NodePositionListInitial,NodesinVeinFramework,TageCVList, Loops)
    for i=1:numel(TageCVList)
        TageCV=TageCVList(i);
        [NumVcp]=FindTargetNumofVe(Loops, TageCV);
        NodePositionListExtracted=ExtractTargetNodes(NodePositionListInitial,TageCV,NumVcp);
        NodePositionListGeneralInitial=GeneralizeNodeNames(NodePositionListExtracted,TageCV,NumVcp);
        [TempVeinFramework]=ObtainVeinFramework(NodesinVeinFramework,NodePositionListGeneralInitial);
        TempVeinFramework.VeinFrame=SpecializeNames(TempVeinFramework.VeinFrame,NumVcp);
        TempVeinFramework.Nodes=SpecializeNames(TempVeinFramework.Nodes,NumVcp);
      
        if i==1
           VeinFramework=TempVeinFramework;
        else
           VeinFramework=cat(1,VeinFramework,TempVeinFramework);
        end
    end
    VeinFrameworkGeometry = unique(VeinFramework,'stable','rows');

end
%%
function [AreabetChordandArc]=ObtainAreabetChordandArc...
    (NodesinVeinFramework,NodePositionListGeneral)
    [NodeNames, NodeCoordinates]=SeparateName(NodePositionListGeneral);
    [FrameNames, NodesofFrame]=SeparateName(NodesinVeinFramework);
    NumFrames=numel(NodesinVeinFramework(:,1));
    FrameLength=zeros(NumFrames,1);

    for i=1:NumFrames
        N=NodesofFrame(i,:);
        NodePositionBox=NodeCoordinates(contains(NodeNames,N),:);
        FrameLength(i)=norm(NodePositionBox(1,:)-NodePositionBox(2,:));
    end

    Ra=FrameLength(contains(FrameNames,"Ra"));
    Rp=FrameLength(contains(FrameNames,"Rp"));
    R=0.5*(Ra+Rp);

    Nea=NodeCoordinates(contains(NodeNames,"Nea"),:);
    Nep=NodeCoordinates(contains(NodeNames,"Nep"),:);
    
    Lchrd=norm(Nea-Nep);
    costheta=(R^2+R^2-Lchrd^2)/(2*R*R);
    theta=acos(costheta);
    sintheta=sin(theta);
    
    AreaTriangle=0.5*R*R*sintheta;
    AreaPie=(theta/(2*pi))*pi*R^2;

    AreabetChordandArc=AreaPie-AreaTriangle;
end



%%
%Calculate Geometric values-------------------------------------------
%--
function [MembraneGeometryOriginal]=ObtainMemGeoOriginal...
    (VeinGeometryOriginal, LoopsOriginal, MembraneAreaPeriList)
    VeinNames=VeinGeometryOriginal.Vein;
    l=VeinGeometryOriginal.Length;
    d=VeinGeometryOriginal.InnerDiameter;
    ContactArea=0.5*pi*l.*d;

    [LoopNames, LoopVeinComs]=SeparateName(LoopsOriginal);
    [LoopVeinComasNum]=TranslateName2Num(LoopVeinComs,VeinNames);
    
    MembraneArea=MembraneAreaPeriList.Area;
    MembraneNames=MembraneAreaPeriList.LoopName;
    NumLoops=numel(LoopNames);
    MaxNumofLoopComponents=(numel(LoopVeinComasNum(1,:)));

    LoopArea=zeros(NumLoops,1);
    Asupply=zeros(NumLoops,1);
    Perimeter=zeros(NumLoops,1);

    for i=1:NumLoops
        f=MembraneNames==LoopNames(i);
        if any(f)
            LoopArea(i)=MembraneArea(f);
            for j=1:MaxNumofLoopComponents
                if ~isnan(LoopVeinComasNum(i,j))
                    Asupply(i)=Asupply(i)+ContactArea(LoopVeinComasNum(i,j));
                    Perimeter(i)=Perimeter(i)+l(LoopVeinComasNum(i,j));
                end
            end 
        end
    end
    
    Aevapo=LoopArea*2;
    Rse=Aevapo./Asupply;
    RAP=sqrt(LoopArea)./Perimeter;
    
    %create table
    sz = [NumLoops,1];
    varNames = "LoopName";
    varTypes = "string";
    
    MembraneGeometryOriginal= table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);
    MembraneGeometryOriginal.LoopName=LoopNames;
    MembraneGeometryOriginal.Rse=Rse;
    MembraneGeometryOriginal.RAP=RAP;
    MembraneGeometryOriginal.Area=LoopArea;
    MembraneGeometryOriginal.Perimeter=Perimeter;
    MembraneGeometryOriginal.SupplyArea=Asupply;
    MembraneGeometryOriginal.EvaporativeArea=Aevapo;
   
end
%--
%Calculate Geometric values-------------------------------------------


%%
function [MembraneGeometrywithNewArea]=ObtainMemArea...
    (MembraneGeometryOld,NodePositionList, ...
    NodesinVeinFramework,TageCV,NumVcp)
    
    MembraneGeometrywithNewArea=MembraneGeometryOld;

    NodesinLe=["Nea";"Nca";"Na";"Np";"Ncp";"Nep";"Nea"];
    NodesinLb=["Nba";"Nca";"Na";"Np";"Ncp";"Nbp";"Nba"];
    
    [ExtractedNodePositionList]=ExtractTargetNodes...
        (NodePositionList, TageCV, NumVcp);

    [NodePositionListGeneralized]=GeneralizeNodeNames...
        (ExtractedNodePositionList, TageCV, NumVcp);

    [NodeNames, NodeCoordinates]=SeparateName(NodePositionListGeneralized);
   
    [AreabetChordandArc]=ObtainAreabetChordandArc...
    (NodesinVeinFramework,NodePositionListGeneralized);
    
    Se=0;
    Sb=0;
    
    Ne=zeros(2:2);
    Nb=zeros(2:2);
    
    for i=1:numel(NodesinLe)-1
        for j=1:2
            Ne(j,:)=NodeCoordinates(contains(NodeNames,NodesinLe(i+j-1)),:);
            %scan coordinates of the target node and its next node for edge
            %membrane
            %first row: coordinates of the target node, second row: the coordinates of next
            %node
            %e.g.,) i=1, j=1: Nea i=1, j=2: Na.....  i=4,j=1: Nep, i=4, j=2: Nea 
    
            Nb(j,:)=NodeCoordinates(contains(NodeNames,NodesinLb(i+j-1)),:);
            %scan coordinates of the target node and its next node for base
            %membane
        end
        
        Se=Se+(Ne(1,1)-Ne(2,1))*(Ne(1,2)+Ne(2,2));
        Sb=Sb+(Nb(1,1)-Nb(2,1))*(Nb(1,2)+Nb(2,2));

    end
    
    EdgeArea=0.5*(abs(Se))+AreabetChordandArc;
    BaseArea=0.5*abs(Sb);
    
    LoopNameEdge=append("L",num2str(NumVcp),"e");
    LoopNameBase=append("L",num2str(NumVcp),"b");

    MembraneGeometrywithNewArea.Area...
        (MembraneGeometrywithNewArea.LoopName==LoopNameEdge)=EdgeArea;
    MembraneGeometrywithNewArea.Area...
        (MembraneGeometrywithNewArea.LoopName==LoopNameBase)=BaseArea;
  
end
%%
function [MembraneGeometrywithNewContact]=ObtainMemContact(LoopsNew,MembraneGeometryOld,VeinGeometryNew)
    
    MembraneGeometrywithNewContact=MembraneGeometryOld;

    VeinNames=VeinGeometryNew.Vein;
    l=VeinGeometryNew.Length;
    d=VeinGeometryNew.InnerDiameter;
    ContactArea=0.5*pi*l.*d;

    [LoopNames, LoopVeinComs]=SeparateName(LoopsNew);
    [LoopVeinComasNum]=TranslateName2Num(LoopVeinComs,VeinNames);

    NumLoops=numel(LoopNames);
    MaxNumofLoopComponents=(numel(LoopVeinComasNum(1,:)));
  
    Asupply=zeros(NumLoops,1);
    Perimeter=zeros(NumLoops,1);

    for i=1:NumLoops
        for j=1:MaxNumofLoopComponents
            if ~isnan(LoopVeinComasNum(i,j))
                Asupply(i)=Asupply(i)+ContactArea(LoopVeinComasNum(i,j));
                Perimeter(i)=Perimeter(i)+l(LoopVeinComasNum(i,j));
            end
        end 
    end

    for scan=1:NumLoops 
        tagerow=MembraneGeometrywithNewContact.LoopName==LoopNames(scan);
        MembraneGeometrywithNewContact.Perimeter(tagerow)=Perimeter(scan);
        MembraneGeometrywithNewContact.SupplyArea(tagerow)=Asupply(scan);

    end
    
end


%%
function [MembraneGeometryNew]=ObtainMembraneGeometry...
    (MembraneGeometryOld,VeinGeometryNew,LoopsNew,NodePositionListNew, NodesinVeinFramework,TageCVList)
    
        
    VeinNames=VeinGeometryNew.Vein;

    Rseo=MembraneGeometryOld.Rse;
    RAPo=MembraneGeometryOld.RAP;
    Amo=MembraneGeometryOld.Area;
    Perio=MembraneGeometryOld.Perimeter;

    MembraneGeometryNew=MembraneGeometryOld;
    for tagecv=1:numel(TageCVList)
        TageCV=TageCVList(tagecv);
        [NumVcp]=FindTargetNumofVe(LoopsNew, TageCV);
        [MembraneGeometryNew]=ObtainMemArea...
            (MembraneGeometryNew,...
            NodePositionListNew, NodesinVeinFramework,...
            TageCV,NumVcp);

    end

    [MembraneGeometryNew]=ObtainMemContact...
            (LoopsNew,MembraneGeometryNew,VeinGeometryNew);
    

    Am=MembraneGeometryNew.Area;
    Peri=MembraneGeometryNew.Perimeter;
    Asupply=MembraneGeometryNew.SupplyArea;

    Aevapo=Am*2;
    Rse=Aevapo./Asupply;
    RAP=sqrt(Am)./Peri;
    
    Rse_Rseo=Rse./Rseo;
    RAP_RAPo=RAP./RAPo;
    Am_Amo=Am./Amo;
    Peri_Perio=Peri./Perio;
    
    MembraneGeometryNew.Rse=Rse;
    MembraneGeometryNew.RAP=RAP;
    MembraneGeometryNew.EvaporativeArea=Aevapo;
    MembraneGeometryNew.Rse_Rseo=Rse_Rseo;
    MembraneGeometryNew.RAP_RAPo=RAP_RAPo;
    MembraneGeometryNew.Am_Amo=Am_Amo;
    MembraneGeometryNew.Peri_Perio=Peri_Perio;

end


%%
function [Nint]=ObtainIntersection(Veins,NodeofVeinList, NodePositionList)
    NodeNames=NodePositionList.Node;
    NodeCoordinates=NodePositionList.Coordinates;
    
    [VeinNames,~]=SeparateName(NodeofVeinList);
    NodesofVein=NodeofVeinList.Nodes;
    
    NodePositionBox=zeros(2,2,2);
    
    for i=1:2
        N=NodesofVein(contains(VeinNames,Veins(i)),:);
        coordinates=NodeCoordinates(contains(NodeNames,N),:);
        NodePositionBox(:,:,i)=coordinates;
    end
    
    system=zeros(2,2);
    solution=zeros(2,1);
    
    for i=1:2
        syspart=(NodePositionBox(1,:,i)-NodePositionBox(2,:,i)).';
        system(:,i)=syspart;
        solpart=(NodePositionBox(2,i,2)-NodePositionBox(2,i,1));
        solution(i)=solpart;
    end
    
    sandt=linsolve(system,solution);
    s=sandt(1);
    %t=sandt(2);
    Nint=s*(NodePositionBox(1,:,1)-NodePositionBox(2,:,1))+NodePositionBox(2,:,1);
end

%%
function [FlowCalcResult]=FlowRateCalculation...
    (VeinGeometry,Loops,FlowDirectionsinLoop,FluidicParameters)
    FlowCalcResult=VeinGeometry;
    VeinNames=FlowCalcResult.Vein;
    
    %--
    %Define fluid parameters-------------------------------------------
    Qin=FluidicParameters.InflowRate;%abs inflow rate
   
    rho=FluidicParameters.Density;%(kg/um^3)Density of hemolymph
    mu=FluidicParameters.Viscosity;%(kg/(um*s))viscousity of hemolymph(Pa*s=kg*m/s^2/m^2*s=kg/(m*s)=10^-6 kg/(um*s))
    g=FluidicParameters.GravitationalAcc;%(um/s^2)gravitational acceleration
    %Define fluid parameters-------------------------------------------
    %--
    
    %--
    %Define vein parameters-------------------------------------------
    d=VeinGeometry.InnerDiameter;%(um)Inner diameters of the wing veins
    l=VeinGeometry.Length;%(um)Lengths of the wing veins
    
    k=128*mu*l./(rho*g*d.^4*pi);%k=h/Q(s/um^2)
    R=128*mu*l./(d.^4*pi);%R=p/Q=(rho*g*h)/Q(s/um^2)
    
    FlowCalcResult.FluidicResistenceHead=k;
    FlowCalcResult.FluidicResistence=R;
    %Define vein parameters-------------------------------------------
    %--
    
    [~,VeinsinLoop]=SeparateName(Loops);
    [~,FlowDirections]=SeparateName(FlowDirectionsinLoop);
    [VeinsinLoopasNums]=TranslateName2Num(VeinsinLoop,VeinNames);
    NumofLoops=numel(VeinsinLoop(:,1));
    LoopNums=1:NumofLoops;
    
    NumofVeins=numel(VeinNames);
    %create a matrix RQ=P system
    Rsystem=zeros(NumofLoops);
    for tageloop=1:NumofLoops
        tagevnum=sum(~isnan(VeinsinLoopasNums(tageloop,:)),2);
        for tagev=1:tagevnum
            TageVeinNum=VeinsinLoopasNums(tageloop,tagev);
            TageLoopNum=LoopNums(logical(sum((VeinsinLoopasNums==TageVeinNum),2)));
            
            Rsystem(tageloop,tageloop)=Rsystem(tageloop,tageloop)+R(TageVeinNum,1);
            for i=TageLoopNum
                if tageloop~=i
                    Rsystem(tageloop,i)=Rsystem(tageloop,i)-R(TageVeinNum,1);           
                end
            end  
        
        end  
    end
    
    Product=Rsystem(:,NumofLoops)*-Qin;
    RsystemNew=Rsystem(:,1:NumofLoops-1);
    Psystem=zeros(NumofLoops,1);
    Psystem(NumofLoops)=-1;
    
    RsystemNew=cat(2,RsystemNew,Psystem);
    Q_P_solutions=linsolve(RsystemNew,Product);
    
    FlowRateforLoops=Q_P_solutions;
    FlowRateforLoops(numel(FlowRateforLoops))=Qin;
    FlowRate=zeros(NumofVeins,1);
    
    for tagevnum=1:NumofVeins
        ContainTageVein=VeinsinLoopasNums==tagevnum;
        ContainTageLoop=logical(sum(ContainTageVein,2));
        if any(ContainTageVein,"all")
            for scan=LoopNums(ContainTageLoop)
                ExtractedDir=FlowDirections(scan,ContainTageVein(scan,:));
                ExtractedFlowrate=FlowRateforLoops(scan);
                FlowRate(tagevnum)=FlowRate(tagevnum)+(ExtractedDir*ExtractedFlowrate);
    
            end
            
        end
    
    end
    
    FlowCalcResult.FlowRate=FlowRate;
    FlowCalcResult.NormalizedFlowRate=FlowRate/Qin;
    FlowCalcResult.PressureLoss_Pa=FlowRate.*R*1e6;
    TotalPressureLoss=nan(NumofVeins,1);
    p=Q_P_solutions(end);
    p_Pa=p*1e6;
    TotalPressureLoss(1)=p_Pa;
    FlowCalcResult.TotalPressureLoss_Pa=TotalPressureLoss;
end