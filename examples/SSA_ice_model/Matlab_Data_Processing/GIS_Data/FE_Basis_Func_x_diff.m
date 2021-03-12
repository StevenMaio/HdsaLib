function [val] = FE_Basis_Func_x_diff(x,y,N,i_index,j_index)
    h = 1/N;
    xi = i_index/N;
    yj = j_index/N;
    val = 0;
    if abs(x-xi)<h && abs(y-yj)<h
        if x < xi
            valx = 1/h;
        else
           valx =  -1/h;
        end
        valy = 1 - (1/h)*abs(y-yj);
        val = valx*valy;
    end
end
