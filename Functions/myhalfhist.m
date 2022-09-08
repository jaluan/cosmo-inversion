function myhalfhist(x0,dx,ymin,ymax,ny,vals,col,lr)

yi = linspace(ymin,ymax,ny); dy=(ymax-ymin)/ny;
[N,~]=histcounts(vals,yi); %[N2,~]=histcounts(uval2,yi);[N3,~]=histcounts(uval3,yi);
Nmax=max(N);
N=dx*N/Nmax; %N2=N2/Nmax; N3=N3/Nmax;
BinCenters=yi(2:end)-0.5*dy;
%patch([smooth(N2(N2>0));zeros(length(N2(N2>0)),1)],[BinCenters(N2>0) fliplr(BinCenters(N2>0))], .7*[1 .1 .1]) %smooth histogram patch
%patch([smooth(N3(N3>0));zeros(length(N3(N3>0)),1)],[BinCenters(N3>0) fliplr(BinCenters(N3>0))], .7*[.3 1 .3], 'Linewidth', 2) %smooth histogram patch

if lr==1
    patch([x0-smooth(N(N>0));x0*ones(length(N(N>0)),1)],[BinCenters(N>0) fliplr(BinCenters(N>0))], col,'EdgeColor','none') %smooth histogram patch
%    patch([x0-f(:);x0*ones(length(f),1)],[yi(:);flipud(yi(:))],col,'linestyle','none');
elseif lr == 2
    patch([x0*ones(length(N(N>0)),1); x0+flipud(smooth(N(N>0)))],[BinCenters(N>0) fliplr(BinCenters(N>0))], col,'EdgeColor','none') %smooth histogram patch
%     patch([x0*ones(length(f),1);x0+flipud(f(:))],[yi(:);flipud(yi(:))],col,'linestyle','none');
end


