function [x1,f1,nf]=powellmain(func,x0,step,Dmat0,tol,maxIt,flag_plot,maxEval)
%[x1,f1,nf]=powellmain(func,x0,step,Dmat0,tol,maxIt,flag_plot,maxEval)
%Powell's method for minimization
%use line scan
%Input:
%   func, function handle
%   x0, initial solution
%   step, step size
%   Dmat0, initial direction set, default to unit vectors
%   tol, a small number to define the termination condition, set to 0 to
%        disable the feature. 
%   maxIt, maximum number of iteration, default to 100
%   flag_plot, 'plot' to plot, otherwise no plot, default to 'noplot'
%   maxEval, maximum number of function evaluation, default to 1500
%Output:
%   x1, best solution
%   f1, func(x1)
%   nf, number of evaluations
%
%Created by X. Huang, 2/22/2013
%[x1,f1,nf]=powellmain(@func_obj,x0,step,dmat)
%[x1,f1,nf]=powellmain(@func_obj,x0,step,dmat,[],100,'noplot')
%[x1,f1,nf]=powellmain(@func_obj,x0,step,dmat,0,100,'noplot',2000)
%
%Reference: X. Huang, et al. Nucl. Instr. Methods, A, 726 (2013) 77-83.
%
%Disclaimer: The RCDS algorithm or the Matlab RCDS code come with absolutely 
%NO warranty. The author of the RCDS method and the Matlab RCDS code does not 
%take any responsibility for any damage to equipments or personnel injury 
%that may result from the use of the algorithm or the code.
%

if nargin <= 2 | isempty(x0)
    step = 0.02;
end
if nargin <= 3 | isempty(Dmat0)
    Dmat0 = eye(length(x0));
end
if nargin <= 4 | isempty(tol)
    tol = 1E-5;
end
if nargin <= 5 | isempty(maxIt)
    maxIt = 100;
end
if nargin <= 6 | isempty(flag_plot)
    flag_plot = 'noplot';
end
if nargin <= 7 | isempty(maxEval)
    maxEval = 1500;
end
Nvar = length(x0);

f0=func(x0);
nf=1;

%best solution so far
xm=x0;
fm=f0;

it=0;
Dmat = Dmat0;

global g_cnt
Npmin = 6;

while it < maxIt
    it = it+1;
    step = step/1.2;
    
    k=1;
    del=0;
    for ii=1:Nvar
        dv = Dmat(:,ii);
        [x1,f1,a1,a2,xflist,ndf] = bracketmin(func,xm,fm,dv,step,flag_plot);
        nf=nf+ndf;
        
        fprintf('iter %d, dir %d: begin\t%d, ',it,ii, g_cnt)
        [x1,f1,ndf] = linescan(func,x1,f1,dv,a1,a2,Npmin,xflist,flag_plot);
        fprintf('end\t%d : %f\n', g_cnt,f1)
        nf=nf+ndf;
        
        %direction with largest decrease
        if (fm - f1)*(1) > del,
            del=(fm - f1)*(1);
            k=ii;
            fprintf('iteration %d, var %d: del = %f updated\n',it, ii, del);
        end
        
        fm=f1;
        xm=x1;
    end
    
    xt = 2*xm-x0;
    ft = func(xt);
    nf = nf+1;
    
    if f0<=ft || 2*(f0-2*fm+ft)*((f0-fm-del)/(ft-f0))^2 >= del %|| norm(xm-x0)<0.01/sqrt(Nvar)
        fprintf('   , dir %d not replaced: %d, %d\n',k,f0<=ft, 2*(f0-2*fm+ft)*((f0-fm-del)/(ft-f0))^2 >= del )
    else
        ndv = (xm-x0)/norm(xm-x0);
        for jj=1:Nvar
            dotp(jj) = abs(ndv(:)'*Dmat(:,jj));
        end
        if max(dotp)<0.9
            if k < Nvar
                Dmat(:,k:Nvar-1)=Dmat(:,k+1:end);
            end
            Dmat(:,end)= ndv;
            
            %move to the mininum of the new direction
            dv = Dmat(:,end);
            [x1,f1,a1,a2,xflist,ndf] = bracketmin(func,xm,fm,dv,step,flag_plot);
            nf=nf+ndf;
            
            fprintf('iter %d, new dir %d: begin\t%d, ',it,k, g_cnt)
            [x1,f1,ndf] = linescan(func,x1,f1,dv,a1,a2,Npmin,xflist,flag_plot);
            fprintf('end\t%d : %f\n',g_cnt,f1)
            nf=nf+ndf;
            fm=f1;
            xm=x1;
        else
            fprintf('    , skipped new direction %d, max dot product %f\n',k, max(dotp));
            
            
        end
        
    end
    
    %termination
    if g_cnt>maxEval
        fprintf('terminated, reaching function evaluation limit: %d > %d\n',g_cnt, maxEval);
        break
    end
    if 2.0*abs(f0-fm) < tol*(abs(f0)+abs(fm)) & tol>0
        fprintf('terminated: f0=%4.2e\t, fm=%4.2e, f0-fm=%4.2e\n',f0, fm, f0-fm);
        break;
    end
    
    f0=fm;
    x0=xm;
end

