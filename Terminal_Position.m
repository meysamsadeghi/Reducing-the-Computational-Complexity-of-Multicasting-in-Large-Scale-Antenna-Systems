function [Terminal_pos] = Terminal_Position (K,G,r) % the inputs are
% r=cell radius and K=number of users
Terminal_pos = zeros(K,G,2);
x = zeros(K,1);
y = zeros(K,1);
%%
for g=1:G % Groups
    for k=1:K
        radious=r*sqrt(rand);
        theta=2*pi*rand;
        x(k,1) = radious*cos(theta) ;
        y(k,1) = radious*sin(theta) ;
        if sqrt(x(k,1)^2 + y(k,1)^2) <= (0.1*r)
            x(k,1) = x(k,1) + (sign(-0.5+rand) * 0.6*r);
            y(k,1) = y(k,1) + (sign(-0.5+rand) * 0.6*r);
        end
    end
    plot(x(:,1),y(:,1),'*', 'Color', [rand,rand,rand]);
    hold on
    Terminal_pos(:,g,1)= x(:,1);
    Terminal_pos(:,g,2)= y(:,1);
end
%
aset = 0:0.1:2*pi;
cell_x = r*sin(aset);
cell_y = r*cos(aset);
hold on
plot(cell_x,cell_y);
hold on
plot(0,0,'s')











