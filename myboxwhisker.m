function myboxwhisker(x0,dx,prc,col)

%whiskers
line([x0,x0],[prc(1),prc(5)],'linewidth',1,'color','k');
line(x0+[-dx,dx],[prc(1),prc(1)],'linewidth',1,'color','k');
line(x0+[-dx,dx],[prc(5),prc(5)],'linewidth',1,'color','k');

%box
patch(x0+[-dx,-dx,dx,dx],[prc(2),prc(4),prc(4),prc(2)],col);

%median
line(x0+[-dx,dx],[prc(3),prc(3)],'linewidth',1,'color','k');
