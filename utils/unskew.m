function om = unskew(S)
% computes the skew-symmetric matrix from the vector om

    if size(S,1) ~= size(S,2)
        error('Invalid skew-symmetric matrix');
    elseif isnumeric(S) && trace(abs(S)) > 1e-10
        error('Invalid skew-symmetric matrix');
    end

    if size(S,1) == 2
        om = S(2,1);
    elseif size(S,1) == 3
        om = zeros(3,1);
        om(1) = S(3,2);
        om(2) = S(1,3);
        om(3) = S(2,1);
    else
        error('bad omega length');
    end
end
