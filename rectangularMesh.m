function [nodeCoordinates, elementNodes] = rectangularMesh(Lx, Ly, numberElementX, numberElementY, elemType)
%RECTANGULARMESH  Generates nodal coordinates and element connectivities for a
%   structured rectangular mesh, considering number of mid-side nodes for Q8 and center node for Q9.
%
%   [nodeCoordinates, elementNodes] = rectangularMesh(Lx, Ly, numberElementX, numberElementY, elemType)
%
%   Inputs:
%       Lx: Length of the rectangle in the x-direction.
%       Ly: Length of the rectangle in the y-direction.
%       numberElementX: Number of elements along the x-direction.
%       numberElementY: Number of elements along the y-direction.
%       elemType:  Element type ('Q4', 'Q8', or 'Q9').
%
%   Outputs:
%       nodeCoordinates: A matrix where each row represents a node, and
%                        the columns are x and y coordinates.
%       elementNodes: A matrix where each row represents an element, and
%                     the columns are the node numbers connected to that
%                     element.  Nodes are numbered going left to right and
%                     bottom to top.

% Validate input arguments
validatestring(elemType, {'Q4', 'Q8', 'Q9'}, 'rectangularMesh', 'elemType');

% Calculate number of nodes based on element type
switch elemType
    case 'Q4'
        numberNodesX = numberElementX + 1;
        numberNodesY = numberElementY + 1;
    case 'Q8'
        numberNodesX = numberElementX * 2 + 1;
        numberNodesY = numberElementY * 2 + 1;
    case 'Q9'
        numberNodesX = numberElementX * 2 + 1;
        numberNodesY = numberElementY * 2 + 1;
        
    otherwise % Add error handling for invalid elemType
      error('Invalid elemType specified. Must be Q4, Q8, or Q9.')
end


% Generate nodal coordinates

switch elemType
    case 'Q4'
        xCoordinates = linspace(0, Lx, numberNodesX);
        yCoordinates = linspace(0, Ly, numberNodesY);
        [X, Y] = meshgrid(xCoordinates, yCoordinates);
        nodeCoordinates = [X(:), Y(:)];
    case 'Q8'
        xCoordinates=linspace(0,Lx,numberNodesX);
        yCoordinates=linspace(0,Ly,numberNodesY);
        [X, Y] = meshgrid(xCoordinates, yCoordinates);
        nodeCoordinates = [X(:), Y(:)];
        
        %some nodes dose not exist and must be deleted
        AY=(2*numberElementY)+1;
        delet_rows=zeros(numberElementX*numberElementY,1);
        AA=zeros(numberElementX,1);
        for j=1:numberElementY
        x=AY+(j*2);
           for i=0:numberElementX-1
           y=x+(2*AY*i);
           AA(i+1,1)=y;
           end
           for k=1:numberElementX
           delet_rows((numberElementY*(k-1))+j,1)=AA(k);
           end
        end
        nodeCoordinates(delet_rows,:)=[];
        
    case 'Q9'
        xCoordinates = linspace(0, Lx, numberNodesX);
        yCoordinates = linspace(0, Ly, numberNodesY);
        [X, Y] = meshgrid(xCoordinates, yCoordinates);
        nodeCoordinates = [X(:), Y(:)];
        
    otherwise % Add error handling for invalid elemType
      error('Invalid elemType specified. Must be Q4, Q8, or Q9.')
end



% Function to calculate the local node locations used in
%   the elementNodes matrix
% Input:
%   j - element row number
%   i - element column number
%   numberNodesX - Number of total nodes in a column
%   elType - Q4, Q8 or Q9 element identifier
% Output:
%   Local node location as calculated according to the Q4, Q8 or Q9.
elementNodes = zeros(numberElementX * numberElementY, getNumNodes(elemType));
elementCounter = 0;

for j = 1:numberElementX
    for i = 1:numberElementY
        elementCounter = elementCounter + 1;

        % Node numbering for Q4 elements
        if strcmp(elemType, 'Q4')
            node1 = (j - 1) * (numberElementY+1) + i;
            node2 = node1 + 1;
            node3 = j * (numberElementY+1) + i + 1;
            node4 = j * (numberElementY+1) + i;
            elementNodes(elementCounter,:) = [node1 node2 node3 node4];
        end

        % Node numbering for Q8 elements
        if strcmp(elemType, 'Q8')

            

                node1 = (j-1)*(3*(numberElementY+1)-1) + 2*(i-1)+1;
                node2 = node1 + 1;
                node3 = node1+2;
                node4 = (j*((numberElementY*2)+1))+j+(j-1)*(numberElementY)+ i;
                node5= (j*(3*(numberElementY+1)-1) + 2*(i-1)+1)+2;
                node6 = node5-1;
                node7 = node5-2;
                node8= node4-1;

          elementNodes(elementCounter,:) = [node1 node2 node3 node4 node5 node6 node7 node8];
        end

        % Node numbering for Q9 elements
        if strcmp(elemType, 'Q9')
             
                node1 = 2*((j-1)*((numberElementY*2)+1))+ 2*(i-1)+1;
                node2 = node1+1;
                node3 = node1+2;
                node4 = (j*((numberElementY*2)+1))+(2*j)+(j-1)*((2*numberElementY)-1)+ (2*i)-1;
                node5= 2*((j)*((numberElementY*2)+1))+ 2*(i-1)+3;
                node6 = node5-1;
                node7= node5-2;
                node8= node4-2;
                node9 =node4-1;
            elementNodes(elementCounter,:) = [node1 node2 node3 node4 node5 node6 node7 node8 node9];


        end

    end
end


function numNodes = getNumNodes(elemType)
 switch elemType
  case 'Q4'
   numNodes = 4;
  case 'Q8'
   numNodes = 8;
  case 'Q9'
   numNodes = 9;
 end
end

end