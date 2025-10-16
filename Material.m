function [AMatrix,BMatrix,DMatrix,SMatrix,Q] = MaterialNew(Ex,Ey,Es,NOOxy,NOOyx,Gxz,Gyz,Teta,Z1,Z2,Kapa)
%This function will calculate the above matrix and Q(local stiffness matrix
%of each layer)

%clacculating Q and QS(stifness matrix for transver strain)
Q11=Ex/(1-NOOxy*NOOyx);
Q22=Ey/(1-NOOxy*NOOyx);
Q12=(NOOxy*Ey)/(1-NOOxy*NOOyx);
Q66=Es;
Q44=Gyz;
Q55=Gxz;
Q=[Q11 Q12 0;Q12 Q22 0;0 0 Q66];
QS=[Q44 0;0 Q55];


%Calculation of Qbar andQbarS(general stiffness matrix of each play)
%and calculating parameters for damping (Qhat and QBighat)

Qbar=zeros(3*size(Teta,2),3);
QbarS=zeros(2*size(Teta,2),2);

for i=1:size(Teta,2)
    m=cosd(Teta(i));
    n=sind(Teta(i));
    
    T=[m^2 n^2 -2*m*n;n^2 m^2 2*m*n;m*n -m*n m^2-n^2];
    TS=[m n;-n m];
    Tinv=inv(T);
    TSinv=inv(TS);

    Qbarloop=Tinv*Q*T;
    QbarSloop=TSinv*QS*TS;

    Qbar([3*i-2:3*i],[1:3])=Qbarloop([1:3],[1:3]);
    QbarS([2*i-1:2*i],[1:2])=QbarSloop([1:2],[1:2]);

end

%calculation of ABD matrix
%and calculating parameters for damping (Ahat,... and ABighat,....)
for i=1:3
    for j=1:3
        A=zeros(3);
        B=zeros(3);
        D=zeros(3);

        for k=1:size(Teta,2)
            a([1:3],[1:3])=Qbar([3*k-2:3*k],[1:3]);

            z1=Z1(1,k);
            z2=Z2(1,k);

            A(i,j)=A(i,j)+(a(i,j)*(z1-z2));
            B(i,j)=B(i,j)+(a(i,j)*(z1^2-z2^2));
            D(i,j)=D(i,j)+(a(i,j)*(z1^3-z2^3));

        end
        AMatrix(i,j)=A(i,j);
        BMatrix(i,j)=(1/2)*B(i,j);
        DMatrix(i,j)=(1/3)*D(i,j);

    end
end



%calculating SMatrix
%and ShatMantrix and SBighatMatrix

SMatrix=zeros(2);
s11=SMatrix(1,1);
s22=SMatrix(2,2);

for k=1:size(Teta,2)
    as([1:2],[1:2])=QbarS([2*k-1:2*k],[1:2]);

    q11=Kapa*as(1,1);
    q22=Kapa*as(2,2);

    z1=Z1(1,k);
    z2=Z2(1,k);
    m=cosd(Teta(k));
    n=sind(Teta(k));


    s11=s11+((q11*m^2)+(q22*n^2))*(z1-z2);
    s22=s22+((q11*n^2)+(q22*m^2))*(z1-z2);

end
SMatrix=[s22 0;0 s11];

end