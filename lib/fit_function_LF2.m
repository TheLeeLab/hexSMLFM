function [model,stdx,mse] = fit_function_LF2(points,terms,pix_per_lens,NA,n,u_scaling,dx)
% FIT_FUNCTION_LF2 fit different lightfield models to list of points passed in array row format (x,y,u,v).
% points:  array row format (x,y,u,v).
% terms: 0 - symettric in u and v, just defocus
%        1 - non-symmetric in u and v, just defocus
%        2 - symettric in u and v, defocus and spherical aberration
%        3 - non symmetric in u and v with defocus and spherical aberration
%        4 - defocus and asymmetric coma
%        5 - defocus and spherical and asymettric coma
%        6 - defocus, spehrical and coma all asymettric
%        7 - an attempt to solve for defocus+spherical in one term to give
%        z*(n/NA)^2
%        8 - analytical solution to the derivative of a spherical wavefront
%        sampled at the centre of a lens position
%        9 - analytical solution averaged across the microlens area
%        10 - allowing u and v directions to have a different alpha value -
%        results in two z values.

u =points(:,2);
v = points(:,3);
x = points(:,4)-pix_per_lens*dx*u;
y = points(:,5)-pix_per_lens*dx*v;
u = u*u_scaling;
v = v*u_scaling;
su = (u.^3+v.^2.*u);
sv = (v.^3+u.^2.*v);
cu = u.^2+(1/3)*v.^2;
cv = v.^2+(1/3)*u.^2;

b = [x;y];
zeros_t = zeros(size(u,1),1);
ones_t = ones(size(u,1),1);

rhosq = (u.^2+v.^2);

switch(terms)
    case 0
        A=[[ones_t zeros_t u]; [zeros_t ones_t v]];
    case 1
        A =[[ones_t zeros_t u zeros_t]; [zeros_t ones_t zeros_t v]];
    case 2
        A = [[ones_t zeros_t u su];[zeros_t ones_t v sv]];
    case 3
        A = [[ones_t zeros_t u zeros_t su zeros_t]; [zeros_t ones_t zeros_t v zeros_t sv]];
    case 4
        A = [[ones_t zeros_t u cu zeros_t]; [zeros_t ones_t v zeros_t cv]];
    case 5
        A =[[ones_t zeros_t u cu zeros_t  su]; [zeros_t ones_t v zeros_t cv sv]];
    case 6
        A =[[ones_t zeros_t u zeros_t cu zeros_t su zeros_t]; [zeros_t ones_t zeros_t v zeros_t cv zeros_t sv]];
    case 7
        A =[[ones_t zeros_t (-(NA/n)^2*u-0.5*(NA/n)^4*su)];[zeros_t ones_t (-(NA/n)^2*v-0.5*(NA/n)^4*sv)]];
    case 8
        A =[[ones_t zeros_t -NA.*u./(n*sqrt(1-rhosq.*NA^2/n^2))];[zeros_t ones_t -NA.*v./(n*sqrt(1-rhosq.*NA^2/n^2))]];
    case 9
         [phi_u, phi_v] = phase_average_sphere(u,v,u_scaling,NA,n,21);
       %[phi_u, phi_v] = phase_average_sphere_integral(u,v,u_scaling,NA,n);
        A =[[ones_t zeros_t phi_u];[zeros_t ones_t phi_v]];   
end


[model,stdx,mse] = lscov(A,b);
