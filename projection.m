clear;

%George Roberts - 03/08/2016 - projection test

u = [2 3 4];
u = u/norm(u);

v = [1 0 0];
w = [0 0 1];

u1 = dot(u,v);
u2 = dot(u,w);

uproj = u1*v + u2*w;
uproj = uproj/norm(uproj);

quiver3(0,0,0,u(1), u(2), u(3),'rx','LineWidth', 2);
hold on;
quiver3(0,0,0,uproj(1), uproj(2), uproj(3),'bx','LineWidth', 2);