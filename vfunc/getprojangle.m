function [ theta_t ] = getprojangle( Q, R0)
%find angle of projection of dipole onto theta-phi plane

R = sqrt(R0(1).^2 + R0(2).^2 + R0(3).^2);
theta = acos(R0(3)./R);
phi = atan2(R0(2),R0(1));

%Get theta unit vector
thx = cos(theta).*cos(phi);
thy = cos(theta).*sin(phi);
thz = -sin(theta);

%Get phi unit vector
phx = -sin(phi);
phy = cos(phi);
phz = 0;

u = dot(Q,[thx thy thz]);
v = dot(Q,[phx phy phz]);

theta_t = atan2d(u,v);

end

