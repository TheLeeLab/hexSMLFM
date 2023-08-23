function dist = forwardModelError(model,points,rhoScaling)

u =points(:,2);
v = points(:,3);
x = points(:,4)-u/rhoScaling;
y = points(:,5)-v/rhoScaling;
alpha = points(:,11:12);

x0 = model(1);
y0 = model(2);
z = model(3);

dist =[(x-x0-z.*alpha(:,1)),(y-y0-z.*alpha(:,2))]; 

end

