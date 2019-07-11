function [x,indxin,indxout]=outlier1d(x,varargin)
%[x,indxin]=outlier1d(x,varargin)
%remove the outliers of x from it
%x is 1d vector with dim>=3
%mul_tol = varargin{1}, default 3.0, 
%perlim  = varargin{2}, default 0.25
%The algorithm is: 1. sort x to ascending order. 2. calculate the difference series of the 
%sorted x. 3. calculate the average difference of the central (1- 2*perlim) part.
%4. examine the difference series of the upper and lower perlim (25% by default) part, if 
%there is a jump that is larger than mul_tol (3 by default) times of the std, remove 
%the upper or lower part from that point on
%
%Created by X. Huang, Nov. 2004
%

mul_tol = 3.0;
len = length(x);
if nargin>=2
	mul_tol = varargin{1}
end
perlim = 0.25;
if nargin>=3
		perlim = varargin{2};
end

if len < 3
	return
elseif len==3 | len==4
	[y,indx]=sort(x);
	dy = diff(y);
	if dy(end) > mul_tol*mean(dy(1:end-1))
		x(indx(end))=[];
	elseif dy(1) > mul_tol*mean(dy(2:end))
		x(indx(1)) = [];
	end
	return
end

%xm = median(x);
[y,indx]=sort(x);
dy = diff(y);

upl = max(floor(len*(1-perlim)), 3);
dnl = max(floor(len*perlim), 2);

%stddy = std(dy(dnl:upl));
stddy = mean(dy(dnl:upl));
upcut = len+1;
dncut = 0;
for ii=upl:len-1
	if dy(ii) > mul_tol * stddy
		upcut = ii+1;
	end
end

for ii=dnl:-1:1
	if dy(ii) > mul_tol * stddy
		dncut = ii;
	end
end
indxin = setxor(1:length(x), indx([1:dncut, upcut:len]));
indxout = indx([1:dncut, upcut:len]);
x(indx([1:dncut, upcut:len])) =[];

