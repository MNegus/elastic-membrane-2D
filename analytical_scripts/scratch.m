x1=linspace(0,20,201);
y1=sin(x1);
x2=linspace(0,20,201);
y2=cos(x2);
sub1=subplot(2,1,1);
axis([0 21 -1.1 1.1]); 
ani1=animatedline('Color','r');
subplot(2,1,2)
axis([0 21 -1.1 1.1]); 
ani2=animatedline('Color','k');
for k=1:length(x1)
    addpoints(ani1,x1(k),y1(k)) ;
    addpoints(ani2,x2(k),y2(k)) ;% I assume you meant to plot (x2,y2)
    drawnow
    pause(0.01);
end