function [Z_grid_out] = Smoothing_Func(Z_grid_in,bandwidth)

    m = size(Z_grid_in,1); % Assumes that Z_grid_in is m x m matrix
    s = zeros(m+2*bandwidth,m+2*bandwidth);
    P1 = 1:bandwidth;
    P2 = (m+bandwidth+1):(m+2*bandwidth);
    I = (bandwidth+1):(bandwidth+m);
    s(I,I) = Z_grid_in;
    s(P1,I) = ones(bandwidth,1)*Z_grid_in(1,:);
    s(P2,I) = ones(bandwidth,1)*Z_grid_in(m,:);
    s(I,P1) = Z_grid_in(:,1)*ones(1,bandwidth);
    s(I,P2) = Z_grid_in(:,m)*ones(1,bandwidth);
    s(P1,P1) = Z_grid_in(1,1);
    s(P1,P2) = Z_grid_in(1,m);
    s(P2,P1) = Z_grid_in(m,1);
    s(P2,P2) = Z_grid_in(m,m);

    Z_grid_out = 0*Z_grid_in;
    for i = (bandwidth+1):(m+bandwidth)
        for j = (bandwidth+1):(m+bandwidth)
            s_loc = s((i-bandwidth):(i+bandwidth),(j-bandwidth):(j+bandwidth));
            Z_grid_out(i-bandwidth,j-bandwidth) = mean(s_loc(:));
        end
    end

end

