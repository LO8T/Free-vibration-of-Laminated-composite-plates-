clear
clc
%profile on -memory

tic;


%Input Geometry:

 Lx=127;
 Ly=12.7;
 Teta=[0 90 0 90 90 0 90 0];
 thickness_ratio=0.08031;   %ratio of thickness to Ly of the plate. (h/Ly)
 A_b = [];
 A_t = [];
for i=1:length(Teta)
     A_b = [A_b -(length(Teta)/2)+i];
     A_t = [A_t -(length(Teta)/2)+(i-1)];
 end
 Z1=((Ly/length(Teta))*thickness_ratio)*A_b; %Cooerdinates of bottom of each layer
 Z2=((Ly/length(Teta))*thickness_ratio)*A_t; %Coordinate of top of each layer
 h=Ly*thickness_ratio;
 thickness=h;
 I=h^3/12;
 numberElementX=20;
 numberElementY=20;
 elemType='Q4';
 typeBC='cfff';


 %Input Material Properties:

 rho=1.48*10^-9;
 Ex=134*10^3;
 Ey=10.3*10^3;
 Es=5*10^3;
 NOOxy=0.33;
 NOOyx=(Ey*NOOxy)/Ex;   %must be calculated remember to move it
 Gxz=3*10^3;
 Gyz=2.5*10^3;


 %Shear correction factor:
 Kapa=(5/6);


 %how many modes I need
 numberOfModes=6;

 
...........................................................................................
 %Calculation Section:

 %calculating ABDS Matrics
 [AMatrix,BMatrix,DMatrix,SMatrix,Q] = Material(Ex,Ey,Es,NOOxy,NOOyx,Gxz,Gyz,Teta,Z1,Z2,Kapa);


 %Creating element and Mesh
 [nodeCoordinates, elementNodes] = rectangularMesh(Lx, Ly, numberElementX, numberElementY, elemType);
 xx=nodeCoordinates(:,1);
 yy=nodeCoordinates(:,2);
 numberElements=numberElementX*numberElementY;


% Nodes and Degree of freedom
numberNodes=size(xx,1);
GDof=5*numberNodes;  %Global number of dgree of freedom


%calculate stiffness and mass matrix
stiffness=formStiffnessMatrixMindlinlaminated5dofLast(GDof,numberElements, ...
    elementNodes,numberNodes,nodeCoordinates,AMatrix,BMatrix,DMatrix,SMatrix,elemType,'complete','reduced');

mass= formMassMatrixMindlinlaminated5dofLast(GDof,numberElements,...
 elementNodes,numberNodes,nodeCoordinates,rho,thickness,Teta,Z1,Z2,I,...
 elemType,'complete');


% boundary conditions
[prescribedDof,activeDof,fixedNodeW] = ...
  EssentialBC5dof( typeBC ,GDof,xx,yy,nodeCoordinates,numberNodes);


%An approach if we get a singular stiffness matrics (adding a fraction of
%mass matrix to both mass and stiffness matrix.
mass=mass+(0*mass);
stiffness=stiffness+(0*mass);


%eigenproblem: free vibrations
% eigenvaluesT is because of manipulation of the mass and stiffness matrices
% in the last section.
[modes,eigenvalues] = eigenvalue(GDof,activeDof,stiffness,mass,numberOfModes);
eigenvaluesT=(eigenvalues*(1+0))-0;
omega=sqrt(eigenvaluesT);
f=omega/(2*pi);


% Liew, p-Ritz 
D0 = Q(2,2)*(h^3/12);  %e2*hˆ3/12/(1-miu12*miu21);
% dimensionless omega
omega_bar = (omega*Ly*Ly/pi^2)*sqrt(rho*h/D0);


% sort out eigenvalues
[omega,ii] = sort(omega);
modes = modes(:,ii);


% drawing mesh and deformed shape
modeNumber = 1;
displacements = modes(:,modeNumber);

w=zeros(numberNodes,1);
for i = 1:numberNodes
    w(i) = displacements(5*(i-1)+3);
end

% surface representation
figure; hold on
for k = 1:size(elementNodes,1)
    patch(nodeCoordinates(elementNodes(k,1:4),1),...
    nodeCoordinates(elementNodes(k,1:4),2),...
    w(elementNodes(k,1:4)),...
    w(elementNodes(k,1:4)))
end
set(gca,'fontsize',18)
view(45,45)


....................................................................................
    %Output Section

display(omega_bar)
display(f)
toc

%profile off;
%p = profile('info') % Get profiler data structure