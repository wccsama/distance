function Xsln = xsolutionCal(BsCoordinate,toaMatrix)
    BsCoordinateSquare    = BsCoordinate.^2;
    toaMatrixSquare = toaMatrix.^2;
    
    [Matrix_m, ~]= size(toaMatrixSquare);
    DimenA = Matrix_m -1;
    A = [];
    pqrParameter =[];
    for i=1:DimenA
        pqrParameter(i)=(BsCoordinateSquare(i,1) + BsCoordinateSquare(i,2)+BsCoordinateSquare(i,3)-toaMatrixSquare(i))...
            -(BsCoordinateSquare(i+1,1) + BsCoordinateSquare(i+1,2)+BsCoordinateSquare(i+1,3)-toaMatrixSquare(i+1));
        Aa = [2*(BsCoordinate(i+1,1)-BsCoordinate(i,1)) 2*(BsCoordinate(i+1,2)-BsCoordinate(i,2)) 2*(BsCoordinate(i+1,3)-BsCoordinate(i,3))];
        A = [A;Aa];
    end
    pqr = -pqrParameter';
    
    Xsln = inv(A'*A)*A'*pqr;
    clear pqrParameter;
    clear pqr;
    clear A;
    clear DimenA;
    clear Matrix_m;
    clear Aa;
end
