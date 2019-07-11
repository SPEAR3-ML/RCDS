function [x1,f1,nf] = linescan(func,x0,f0,dv,alo,ahi,Np,xflist,flag_plot)
%line scan in the parameter space along a direction dv
%Input:
%   func, function handle, f=func(x)
%   x0,  Npx1 vec, initial point
%   f0,  f0=func(x0), provided so that no need to re-evaluate, can be NaN
%        or [], to be evaluated. 
%   dv, Npx1  vec, the unit vector for a direction in the parameter space
%   alo, ahi, scalar, the low and high bound of alpha, for f=func(x0+a dv)
%   Np, minimum number of points for fitting. 
%   xflist, Nx2, known solutions
%Output:
%   x1, the new solution
%   f1, the value at the new solution, f1=func(x1)
%   nf, the number of function evaluations
%Created by X. Huang, 1/25/2013
%
if nargin<7
    flag_plot = 1;
end
nf =0;
if isnan(f0) | isempty(f0)
    f0 = func(x0);
    nf = nf + 1;
end

if alo>=ahi
    error('high bound should be larger than the low bound');
end
V = length(x0);
if length(x0)~=length(dv)
    error('parameter space dimension not consistent');
end
if isempty(Np) 
    Np = 6;
end
if Np<6
    Np=6;
end
delta = (ahi-alo)/(Np-1);


alist = alo:delta:ahi;
a0indx = find(xflist(:,1)>=alo & xflist(:,1)<=ahi);
a0list = xflist(a0indx,1)';
f0list = xflist(a0indx,2)';
indx_rmv = [];
for ii=1:length(alist)
    if min(abs(alist(ii)-a0list))<= delta/2.0
       indx_rmv = [indx_rmv,ii]; 
    end
end
alist(indx_rmv) =[];

if 1 %evaluate here
    flist = zeros(size(alist))*NaN;
    for ii=1:length(alist)
        alpha = alist(ii);
        flist(ii) = func(x0+alpha*dv);
        nf = nf+1;
    end
    
else %use mass_evaluate
    xa = zeros(length(alist), V);
    for ii=1:length(alist)
        alpha = alist(ii);
        xa(ii,:) = x0+alpha*dv;
    end
    
    xa = mass_evaluate(xa,1,V,func);
    flist = xa(:,V+1)';
    nf = nf+V;
end
alist = [alist, a0list];
flist = [flist, f0list];
[alist, sindx] = sort(alist);
flist = flist(sindx);

% filter out NaNs
[tmpNan] = isnan(flist);
alist(tmpNan) = [];
flist(tmpNan) = [];
if length(alist)<=0
    x1=x0;
    f1 = func(x0);
    nf = nf+1;
    return
end

[fn1,imo] = min(flist);
xn1 = x0+alist(imo)*dv;
if length(alist)<=5
    x1 = xn1;
    f1 = fn1;
    return;   
end

%
MP = 101;
tmp_min = max(alist(1), alist(imo)-6*delta);
tmp_max = min(alist(end), alist(imo)+6*delta);
av = tmp_min + (tmp_max-tmp_min)*(0:MP-1)*1.0/MP;

if 1
     warning('off', 'MATLAB:polyfit:RepeatedPointsOrRescale');
    %parabolic interpolation
    [p,s] = polyfit(alist, flist,2);
    cfl = polyval(p,alist);
    [ftmp,indxin,indxout]=outlier1d(flist-cfl);
    if length(indxout)<=1
        if length(indxout)==1
            [p,s] = polyfit(alist(indxin), flist(indxin),2);
        end
        yv = polyval(p, av);
        [myv,iyv] = min(yv);
        
        x1 =x0+av(iyv)*dv;
        % f1 = func(x1);
        f1 = myv;
    
    else %no interpolation
       x1 = xn1;
       f1 = fn1;
       yv = interp1(alist,flist,av);
           [myv,iyv] = min(yv);

    end
    % flag_plot = 1;
    if strcmp(flag_plot,'plot')
        figure(1210)
        plot(alist, flist,'o-',av,yv,'-',0,f0,'r*',av(iyv),f1,'rv');
        set(gca,'xlim',[alo-delta,ahi+delta]);
        xlabel('alpha')
        ylabel('objective')
        global g_cnt
        saveas(1210,['pl_fit_' num2str(g_cnt) '.fig']); 
    end

else
    %interpolate
    % yv = interp1(alist,flist,av,'spline');
    yv = interp1(alist,flist,av);
    
    [myv,iyv] = min(yv);
    x1 =x0+av(iyv)*dv;
    % f1 = func(x1);
    f1 = myv;
end


