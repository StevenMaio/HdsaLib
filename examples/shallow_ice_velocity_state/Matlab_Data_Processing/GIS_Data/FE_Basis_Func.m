function [val] = FE_Basis_Func(x,y,N,i_index,j_index)
    h = 1/N;
    xi = i_index/N;
    yj = j_index/N;
    val = 0;
    if abs(x-xi)<h && abs(y-yj)<h
        %valx = 1 - (2/h^2)*(x-xi)^2 + (1/h^4)*(x-xi)^4;
        %valy = 1 - (2/h^2)*(y-yj)^2 + (1/h^4)*(y-yj)^4;
        valx = 1 - (1/h)*abs(x-xi);
        valy = 1 - (1/h)*abs(y-yj);
        val = valx*valy;
    end
end

