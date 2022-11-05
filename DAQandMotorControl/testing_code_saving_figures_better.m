%%
x = (0 : 0.01 : 1)';
y = (0 : 0.01 : 1)';
[X, Y] = meshgrid(x, y);
Z = peaks(X, Y);
Z(40:60, 60:80) = NaN;
figure
contourf(X, Y, Z)
colorbarpwn(min(Z, [], 'all'), max(Z, [], 'all'))
ax = gca;
ax.Color = 'r';
