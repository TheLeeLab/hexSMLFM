function [alpha_u,alpha_v] = phase_average_sphere_integral(u,v,u_scaling,NA,n)
d=u_scaling;
alpha_u=n/(d*NA)*(sqrt(1-(NA^2/n^2)*((u+d/2).^2+v.^2))-sqrt(1-(NA^2/n^2)*((u-d/2).^2+v.^2)));

alpha_v=n/(d*NA)*(sqrt(1-(NA^2/n^2)*((v+d/2).^2+u.^2))-sqrt(1-(NA^2/n^2)*((v-d/2).^2+u.^2)));

end

