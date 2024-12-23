function [val] = FE_Basis_Func_y_diff(x,y,N,i_index,j_index)
    h = 1/N;
    xi = i_index/N;
    yj = j_index/N;
    val = 0;
    if abs(x-xi)<h && abs(y-yj)<h
        if y < yj
            valy = 1/h;
        else
           valy =  -1/h;
        end
        valx = 1 - (1/h)*abs(x-xi);
        val = valx*valy;
    end
end
