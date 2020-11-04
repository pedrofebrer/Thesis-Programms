function err = RMSE(Res,X)
    
    % Input:
    % err - mean root square error value
    
    % Output:
    % Res - Results two dimensional matrix
    % X - Formula\Simulation matrix values
    
    v = size(Res);  % v is defined as the vector with the width and length of matrix Res as its values
    s = 0;          % s is defined as 0
    n = 0;          % n is defined as 0
    
    for i = 1:v(1)
        for j = 1:v(2)
            if(Res(i,j) ~= 0)
                s = s + (Res(i,j)-X (i,j))^2;   % s is increased by the square value of the difference of the (i,j)-elements of matrices Res and X
                n=n+1;
            end
        end
    end
    
    err = sqrt(s/n);    % err is defined as the square root of of the mean of the square difference of the values of matrices Res and X
                        % i.e. the root mean square error value
end