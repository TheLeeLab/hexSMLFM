function [alpha_u,alpha_v] = phase_average_sphere(u,v,u_scaling,NA,n,m)

du = u_scaling/m;
du2 = u_scaling/2;
alpha_u = u;
alpha_v = v;

for i = 1:numel(u)
    u_min = (u(i)-du2);
    u_max = (u(i)+du2);
    v_min = (v(i)-du2);
    v_max = (v(i)+du2);
    [um, vm] = meshgrid(u_min:du:u_max,v_min:du:v_max);
    rho = sqrt(um.^2+vm.^2);
    
    
    phi_u = -(NA/n).*um./(sqrt(1-(NA/n)^2.*rho.^2));
    phi_u(imag(phi_u)~=0) = NaN;
    alpha_u(i) = nanmean(phi_u(:));
    
    phi_v = -(NA/n).*vm./(sqrt(1-(NA/n)^2.*rho.^2));
    phi_v(imag(phi_v)~=0) = NaN;
    alpha_v(i) = nanmean(phi_v(:));
    
end

end

