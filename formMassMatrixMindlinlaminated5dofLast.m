%In this function The M will be calculated with parameters order of
%[u v w fix fiy]

function [M]= ...
 formMassMatrixMindlinlaminated5dofLast(GDof,numberElements,...
 elementNodes,numberNodes,nodeCoordinates,rho,thickness,Teta,Z1,Z2,I,...
 elemType,quadType)

 % computationof mass matrix
 % for Mindlinplate element

 Ip=rho*thickness;
 Ic=0;
 Ib=0;
 for i=1:length(Teta)
     Ib=Ib+(rho/3)*(Z1(i)^3-Z2(i)^3);
 end

 M1=[Ip 0 0 0 0;0 Ip 0 0 0;0 0 Ip 0 0;0 0 0 Ib 0;0 0 0 0 Ib];
 M2=[0 0 0 -Ic 0;0 0 0 0 -Ic;0 0 0 0 0;-Ic 0 0 0 0;0 -Ic 0 0 0];

 M = zeros(GDof);

 % Gauss quadraturefor bendingpart
 [gaussWeights,gaussLocations]= gaussQuadrature(quadType);

 % cycle forelement
 for e = 1:numberElements
     % indice: nodal connectivities for each element
     % elementDof: element degrees of freedom
     indice = elementNodes(e,:);

     elementDof=zeros(1,5*size(elementNodes,2));
     for i=1:size(elementNodes,2)
         elementDof((5*i)-4:5*i)=[(indice(i)*5)-4:indice(i)*5];
     end
     ndof = length(indice);
     % cycle for Gauss point

     for q = 1:size(gaussWeights,1)
         GaussPoint = gaussLocations(q,:);
         xi = GaussPoint(1);
         eta = GaussPoint(2);


         % shape functions and derivatives
         [shapeFunction,naturalDerivatives] = ...
         shapeFunctionsQ(xi,eta,elemType);
         % Jacobian matrix, inverse of Jacobian,
         % derivatives w.r.t. x,y
         [Jacob,invJacobian,XYderivatives] = ...
         Jacobian(nodeCoordinates(indice,:),naturalDerivatives);


         % [B] matrix shear
         N = zeros(5,5*ndof);
         for i=1:ndof
             N(1,(5*i)-4)=shapeFunction(i);
             N(2,(5*i)-3)=shapeFunction(i);
             N(3,(5*i)-2)=shapeFunction(i);
             N(4,(5*i)-1)=shapeFunction(i);
             N(5,(5*i)-0)=shapeFunction(i);
         end

         % stiffness matrix shear
         M(elementDof,elementDof) = M(elementDof,elementDof) + ...
         N'*(M1+M2)*N*gaussWeights(q)*det(Jacob);

     end % end gauss point loop
 end % end element loop
 end