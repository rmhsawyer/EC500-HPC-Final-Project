x = rand([15 2]);
y = zeros(15, 1);
for i = 1:15
    y(i) = sign([2 -1]* x(i,:)' - 1);
end
lx = 0.5:0.1:1;
ly = 2.*lx-1;

figure
hold on
plot(x(y==1,1), x(y==1,2),'b*')
plot(x(y==-1,1), x(y==-1,2),'r*')
plot(lx,ly)
hold off
