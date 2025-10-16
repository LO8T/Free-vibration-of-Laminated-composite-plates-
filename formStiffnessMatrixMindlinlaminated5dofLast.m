%In this function The K will be calculated with parameters order of
%[u,v,w,fix,fiy]

function [K] = ...
formStiffnessMatrixMindlinlaminated5dofLast(GDof,numberElements,elementNodes,numberNodes,nodecoordinates,AMatrix,BMatrix,DMatrix,SMatrix,elemType,quadTypeB,quadTypeS)
% computation of stiffness matrix for laminated plate element
% K: stiffness matrix

K = zeros(GDof);

% Gauss quadrature for bending part
[gaussWeights,gaussLocations] = gaussQuadrature(quadTypeB);

% cycle for element
for e = 1:numberElements
    % indice: nodal connectivities for each element
    % elementDof: element degrees of freedom
    indice = elementNodes(e,:);
    elementDof=zeros(1,5*size(elementNodes,2));
    for i=1:length(indice)
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
        Jacobian(nodecoordinates(indice,:),naturalDerivatives);


        % [B] matrix bending
        B_f = zeros(3,5*ndof);
        for i=1:ndof
            B_f(1,(5*i)-1)=XYderivatives(i,1);
            B_f(2,(5*i)-0)=XYderivatives(i,2);
            B_f(3,(5*i)-1)=XYderivatives(i,2);
            B_f(3,(5*i)-0)=XYderivatives(i,1);
        end

        % [B] matrix 
        B_p = zeros(3,5*ndof);
        for i=1:ndof
            B_p(1,(5*i)-4)=XYderivatives(i,1);
            B_p(2,(5*i)-3)=XYderivatives(i,2);
            B_p(3,(5*i)-4)=XYderivatives(i,2);
            B_p(3,(5*i)-3)=XYderivatives(i,1);
        end


        % stiffness matrix
        % ... bending-bending
        K(elementDof,elementDof) = K(elementDof,elementDof) + ...
        B_f'*DMatrix*B_f*gaussWeights(q)*det(Jacob);
        % ... membrane-membrane
        K(elementDof,elementDof) = K(elementDof,elementDof) + ...
        B_p'*AMatrix*B_p*gaussWeights(q)*det(Jacob);
        % ... membrane-bending
        K(elementDof,elementDof) = K(elementDof,elementDof) + ...
        B_p'*BMatrix*B_f*gaussWeights(q)*det(Jacob);
        % ... bending-membrane
        K(elementDof,elementDof) = K(elementDof,elementDof) + ...
        B_f'*BMatrix*B_p*gaussWeights(q)*det(Jacob);


    end % Gauss point
end % element



% shear stiffness matrix
% Gauss quadrature for shear part
[gaussWeights,gaussLocations] = gaussQuadrature(quadTypeS);

% cycle for element
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
        Jacobian(nodecoordinates(indice,:),naturalDerivatives);


        % [B] matrix shear
        B_s = zeros(2,5*ndof);
        for i=1:ndof
            B_s(1,(5*i)-2)=XYderivatives(i,1);
            B_s(2,(5*i)-2)=XYderivatives(i,2);
            B_s(1,(5*i)-1)=shapeFunction(i);
            B_s(2,(5*i)-0)=shapeFunction(i);
        end

        % stiffness matrix shear
        K(elementDof,elementDof) = K(elementDof,elementDof) + ...
        B_s'*SMatrix*B_s*gaussWeights(q)*det(Jacob);

    end % end gauss point loop
end % end element loop


end