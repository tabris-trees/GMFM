function kernel_matrix_station = KM(phi,theta,order)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
    syms dirac_0;
    dirac_0 = sym(0);
       
    cos_theta = cos(theta);
    sin_theta = sin(theta); 
    
    kernel_matrix_station = [];
    kernel_M_X = [];
    kernel_M_Y = [];
    kernel_M_Z = [];
    legendre_P_n0 = 1; % P00 = 1
    
    for n = 1:order
        m = 0:n;
        n_mat = ones(length(m),1)'*n;
        two_mat = ones(length(m),1)'*2;
        dirac_m = double(kroneckerDelta(m,dirac_0));
%         Factor_norm = (-1).^m.*sqrt((two_mat-dirac_m).*(factorial(n_mat-m)./factorial(n_mat+m)));
        Factor_norm = sqrt((two_mat-dirac_m).*(factorial(n_mat-m)./factorial(n_mat+m)));
        legendre_P_n = Factor_norm.*legendre(n,cos_theta)';
        D_legendre_P_n = (n*cos_theta.*legendre_P_n-sqrt(n_mat.^2-m.^2).*[legendre_P_n0,0])./sin_theta;
        g_kernel_X = cos(phi*m).*D_legendre_P_n;
        h_kernel_X = sin(phi*m).*D_legendre_P_n;
        kernel_M_X = [kernel_M_X,g_kernel_X,h_kernel_X];
        g_kernel_Y = m.*sin(phi*m).*legendre_P_n./sin_theta;
        h_kernel_Y = -m.*cos(phi*m).*legendre_P_n./sin_theta;
        kernel_M_Y = [kernel_M_Y,g_kernel_Y,h_kernel_Y];
        g_kernel_Z = -cos(phi*m).*legendre_P_n.*(n+1);
        h_kernel_Z = -sin(phi*m).*legendre_P_n.*(n+1);
        kernel_M_Z = [kernel_M_Z,g_kernel_Z,h_kernel_Z];
    
        legendre_P_n0 = legendre_P_n;
    end
    % Remove the zero element (hn0)
    index_remove = (2:order+1).^2-ones(1,order);
    kernel_M_X(index_remove) = [];
    kernel_M_Y(index_remove) = [];
    kernel_M_Z(index_remove) = [];
    
    kernel_matrix_station = [kernel_matrix_station;kernel_M_X;kernel_M_Y;kernel_M_Z];
end