function [dist] = dist_fcn_LF2(model,points,terms,pix_per_lens,NA,n,u_scaling,dx)

% terms: 0 - symettric in u and v, just defocus
%        1 - non-symmetric in u and v, just defocus
%        2 - symettric in u and v, defocus and spherical aberration
%        3 - non symmetric in u and v with defocus and spherical aberration
%        4 - defocus and asymmetric coma
%        5 - defocus and asymmetric coma and symmetric spherical
%        6 - defocus, spehrical and coma all asymettric

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
rhosq = (u.^2+v.^2);
x0 = model(1);
y0 = model(2);


switch(terms)
    case 0
        a = model(3);
        dist = [x-(x0+u.*a),y-(y0+v.*a)];
    case 1
        a = model(3);
        b = model(4);
        dist = [(x-(x0+u.*a)),(y-(y0+v.*b))];
    case 2
        a = model(3);
        c = model(4);
        dist = [(x-(x0+u.*a+c*(u.^3+v.^2.*u))),(y-(y0+v.*a+c*(v.^3+u.^2.*v)))];
    case 3
        a = model(3);
        b = model(4);
        c = model(5);
        d = model(6);
        dist = [x-(x0+u.*a+c*(u.^3+v.^2.*u)),y-(y0+v.*b+d*(v.^3+u.^2.*v))];
    case 4
        a = model(3);
        b = model(4);
        c = model(5);
        dist =  [(x-(x0+u.*a+b*(3*u.^2+v.^2))),+(y-(y0+v.*a+c*(3*v.^2+u.^2)))];
    case 5
        a = model(3);
        b = model(4);
        c = model(5);
        d = model(6);
        dist =  [(x-(x0+u.*a+b*(3*u.^2+v.^2)+d*(u.^3+v.^2.*u))),+(y-(y0+v.*a+b*(3*v.^2+u.^2)+c*(v.^3+u.^2.*v)))];
    case 6
        a = model(3);
        b = model(4);
        c = model(5);
        d = model(6);
        e = model(7);
        f = model(8);
        dist =  [(x-(x0+u.*a+c*(3*u.^2+v.^2)+e*(u.^3+v.^2.*u))),+(y-(y0+v.*b+d*(3*v.^2+u.^2)+f*(v.^3+u.^2.*v)))];
    case 7
        b = model(3);
        dist = [(x-x0-b*(-(NA/n)^2*u-0.5*(NA/n)^4*su)),(y-y0-b*k*(-(NA/n)^2*v-0.5*(NA/n)^4*sv))];
    case 8
        b = model(3);
        alpha_u=-NA*u./(n*sqrt(1-rhosq.*NA^2/n^2));
        alpha_v=-NA*v./(n*sqrt(1-rhosq.*NA^2/n^2));
        dist = [(x-x0-b*alpha_u),(y-y0-b*alpha_v)];
    case 9
        b = model(3);
        [phi_u, phi_v] = phase_average_sphere(u,v,u_scaling,NA,n,21);
        %[phi_u, phi_v] = phase_average_sphere_integral(u,v,u_scaling,NA,n);
        dist =[(x-x0-b*phi_u),(y-y0-b*phi_v)]; 
end

