function [SquareH] = spaRegMatrix(N,spatialRreg)
% SquareH - regularization operation
% spatialRreg -regularization method: 'zeros-spatial-reg' '1st-spatial-reg' '2nd-spatial-reg'
% N - number of Degrees of freedom (DOF)
if N > 3
    SquareH =zeros(N,N);
    switch spatialRreg
        case '2nd-spatial-reg'
            SquareH(1,1)=2;  SquareH(1,2)=-4;  SquareH(1,3)=2;
            SquareH(2,1)=-4; SquareH(2,2)=10;  SquareH(2,3)=-8; SquareH(2,4)=2;
            j=1;
            for i=3:(N-2)
                SquareH(i,j)=2;     SquareH(i,j+1)=-8;
                SquareH(i,j+2)=12;
                SquareH(i,j+3)=-8;  SquareH(i,j+4)=2;
                j=j+1;
            end
            SquareH(N-1,N-3)=2; SquareH(N-1,N-2)=-8;  SquareH(N-1,N-1)=10; SquareH(N-1,N)=-4;
            SquareH(N,N-2)=2;  SquareH(N,N-1)=-4;  SquareH(N,N)=2;
            SquareH =SquareH./2;
        case '1st-spatial-reg'
            SquareH(1,1)=1;SquareH(1,2)=-1;
            j=1;
            for i=2:(N-1)
                SquareH(i,j)= -1;
                SquareH(i,j+1)= 2;
                SquareH(i,j+2)= -1;
                j=j+1;
            end
            SquareH(N,N-1)=-1;SquareH(N,N)=1;
        case 'zeros-spatial-reg'
            for i=1:N
                SquareH(i,i)= 1;
            end
        otherwise
            disp('spatial_reg: input spatial regularization');
    end
else
    disp('error, N - Degrees of freedom must larger than 3')
end
end