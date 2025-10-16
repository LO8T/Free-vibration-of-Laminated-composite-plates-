function [JacobianMatrix,invJacobian,XYDerivatives] = ...
Jacobian(nodeCoordinates,naturalDerivatives)
% JacobianMatrix: Jacobian matrix
% invJacobian: inverse of Jacobian Matrix
% XYDerivatives: derivatives w.r.t. x and y
% naturalDerivatives: derivatives w.r.t. xi and eta
% nodeCoordinates: nodal coordinates at element level
JacobianMatrix = nodeCoordinates'*naturalDerivatives;    %it is in 9.2.48 ready
invJacobian = inv(JacobianMatrix);
XYDerivatives = naturalDerivatives*invJacobian;   % it is in 9.2.40a ready
end % end function Jacobian