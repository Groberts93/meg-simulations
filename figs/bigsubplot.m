% Load saved figures
dr=hgload('dipole_r.fig');
dt=hgload('dipole_t.fig');
dp=hgload('dipole_p.fig');

vr=hgload('volume_r.fig');
vt=hgload('volume_t.fig');
vp=hgload('volume_p.fig');

% Prepare subplots
figure
for i = 1:6
    h(i) = subplot(2,3,i);
    view(3);
end
% Paste figures on the subplots
copyobj(allchild(get(dr,'CurrentAxes')),h(1));
view(3);
copyobj(allchild(get(dt,'CurrentAxes')),h(2));
view(3);
copyobj(allchild(get(dp,'CurrentAxes')),h(3));
view(3);
copyobj(allchild(get(vr,'CurrentAxes')),h(4));
view(3);
copyobj(allchild(get(vt,'CurrentAxes')),h(5));
view(3);
copyobj(allchild(get(vp,'CurrentAxes')),h(6));
view(3);


