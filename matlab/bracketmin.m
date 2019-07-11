function [xm,fm,a1,a2,xflist,nf] = bracketmin(func,x0,f0,dv,step,flag_plot)
%bracket the minimum along the line with unit direction dv
%Input:
%   func, function handle, f=func(x)
%   x0,  Npx1 vec, initial point, func(x0+alpha*dv)
%   f0,  f0=func(x0), provided so that no need to re-evaluate, can be NaN
%        or [], to be evaluated. 
%   dv, Npx1  vec, the unit vector for a direction in the parameter space
%   step, initial stepsize of alpha,
%Output:
%   xm, fm, the best solution and its value
%   a1,a2, values of alpha that satisfy f(xm+a1*dv)>fmin and
%   f(xm+a2*dv)>fmin
%   xflist, Nx2, all tried solutions
%   nf, number of function evluations
%created by X. Huang, 1/25/2013
%

if nargin<5
    flag_plot = 1;
end
nf = 0;

global g_noise
if isempty(g_noise) %isempty(whos('global','g_noise'))
    g_noise = 0.0;
end
if isnan(f0) | isempty(f0)
    f0 = func(x0);
    nf = nf + 1;
end

if 1
    xflist(1,:) = [0,f0];
end
fm = f0;
am = 0;
xm = x0;

step_init = step;

x1 = x0+dv*step;
f1 = func(x1);
nf = nf + 1;
if 1
    nxf = size(xflist,1)+1;
    xflist(nxf,:) = [step,f1];
end
if f1<fm
   fm = f1;
   am = step;
   xm = x1;
end

gold_r = 1.618034;
while f1<fm+g_noise*3
    if abs(step)<0.1  %maximum step
        step = step*(1.0+gold_r);
    else
        step = step+0.1;
    end
    x1 = x0+dv*step;

    f1 = func(x1);
    nf = nf + 1;
    if isnan(f1)
        step = step/(1.0+gold_r); %get the last valid solution
        fprintf('bracketmin: f1=NaN\n')
       break; 
    end
    if 1
        nxf = size(xflist,1)+1;
        xflist(nxf,:) = [step,f1];
    end
    if f1<fm
        fm = f1;
        am = step;
        xm = x1;
    end
   
end
a2 = step;
if f0>fm+g_noise*3
   a1=0;

   if strcmp(flag_plot,'plot')
       figure(1302)
       plot(xflist(:,1), xflist(:,2),'o',xflist(1,1), xflist(1,2),'r*');
       
   end
   a1 = a1-am;
   a2 = a2-am;
   xflist(:,1) = xflist(:,1)-am;
   return;
end

step = -step_init;
x2 = x0+dv*step;
f2 = func(x2);
nf = nf + 1;
if 1
    nxf = size(xflist,1)+1;
    xflist(nxf,:) = [step,f2];
end
if f2<fm
    fm = f2;
    am = step;
    xm = x2;
end

while f2<fm+g_noise*3
    if abs(step)<0.1  %maximum step
        step = step*(1.0+gold_r);
    else
       step = step-0.1; 
    end
    x2 = x0+dv*step;
    f2 = func(x2);
    nf = nf + 1;
    if isnan(f2)
        step = step/(1.0+gold_r);
        fprintf('bracketmin: f2=NaN\n')
        break;
    end
    
    if 1
        nxf = size(xflist,1)+1;
        xflist(nxf,:) = [step,f2];
    end
    if f2<fm
        fm = f2;
        am = step;
        xm = x2;
    end

end
a1 = step;

if a1>a2
    tmpa = a2;
    a2 = a1;
    a1=tmpa;
end
a1 = a1-am;
a2 = a2-am;
xflist(:,1) = xflist(:,1) - am;

[sa,si] = sort(xflist(:,1));
xflist = xflist(si,:);
if strcmp(flag_plot,'plot')
    figure(1302)
    plot(xflist(:,1), xflist(:,2),'o',-am, f0,'r*',0,fm,'rd');
    xlabel('alpha')
    ylabel('objective');
    
end
