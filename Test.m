x = linspace(0, 5);
y = linspace(0, 10);
[x, y] = meshgrid(x, y); 

z = exp(-y.*x).*(y <= 5) + exp((y-10).*x).*(y > 5); 

surf(x, y, z)