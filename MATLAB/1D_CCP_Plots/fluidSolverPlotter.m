fileID = fopen('sin3.txt');
data = textscan(fileID, '%f %f %f %f', 'HeaderLines',1);
x = data{1}; t = data{2}; u_approx = data{3}; u_exact = data{4};

n_cell = 64;
ts = length(x)/n_cell;

x = reshape(x, n_cell, ts);
t = reshape(t, n_cell, ts);
u_approx = reshape(u_approx, n_cell, ts);
u_exact = reshape(u_exact, n_cell, ts);

figure
pos = get(gcf, 'position');
set(gcf, 'Position', [pos(1), pos(2), pos(1) + 3*(pos(3)-pos(1)), pos(4)]);
subplot(1,2,1)
h = surf(x,t,u_approx);
get(h);
set(h,'linestyle','none');
title('Approximation')
xlabel('x')
ylabel('t')
zlabel('u')

subplot(1,2,2)
h = surf(x,t,u_exact);
get(h);
set(h,'linestyle','none');
title('Exact')
xlabel('x')
ylabel('t')
zlabel('u')


