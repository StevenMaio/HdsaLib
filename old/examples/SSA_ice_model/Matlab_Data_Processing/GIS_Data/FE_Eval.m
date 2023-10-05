function [val] = FE_Eval(x,y,c,N)
    i = floor(x*N);
    j = floor(y*N);
    val = 0;
    for ii = -1:1
       for jj = -1:1
            i_index = i+ii;
            j_index = j+jj;
            if i_index>=0 && j_index>=0 && i_index<N+1 && j_index<N+1
               k = j_index*(N+1)+i_index + 1; % Remove +1 in C++ 
               val = val + c(k)*FE_Basis_Func(x,y,N,i_index,j_index);
            end
       end
    end
end

