%Disclaimer: The RCDS algorithm or the Matlab RCDS code come with absolutely 
%NO warranty. The author of the RCDS method and the Matlab RCDS code does not 
%take any responsibility for any damage to equipments or personnel injury 
%that may result from the use of the algorithm or the code.
%

clear

global vrange Nvar

Nobj = 1;
Nvar = 6;

vrange = [-ones(Nvar,1), ones(Nvar,1)]*150; 

%initial solution, in this test it is random
p0 = randn(1,Nvar)*5; 

%in reality, we may want to read the present setpoints as below
%p0 = getSetPt;  %write your function getSetPt

x0 = (p0'-vrange(:,1))./(vrange(:,2)-vrange(:,1));

global g_cnt g_data
g_data=[];
g_cnt = 0;

%% evaluate the random noise level of the objective function
% no need to evaluate every time 
for ii=1:100
   tmp_data(ii) = func_obj(x0); 
end

global g_noise;
%g_noise = 0.001;  
g_noise = std(tmp_data);  

%% 

% dfname = appendtimestamp('diary');
% diary(dfname);

%use unit vectors as initial direction set
dmat = eye(Nvar);

step = 0.01;

[x1,f1,nf]=powellmain(@func_obj,x0,step,dmat,0,100,'noplot',1000)
% diary off

save tmpRCDSdata

% N = size(g_data,1);
% figure
% plot(1:N, g_data(:,Nvar+1),'.')

%% save data
% descrip = 'tmp1';
% eval(['mkdir ' descrip])
% % eval(['!mv pl_fit*.fig ' descrip])
% eval(['!mv diary* ' descrip])
% eval(['!mv tmpRCDSdata.mat  data_' descrip '.mat']);
% eval(['!cp  data_' descrip '.mat  ' descrip]);
% eval(['!cp  *.m  ' descrip]);

[data_1,xm,fm] = process_scandata(g_data,Nvar,vrange,'plot');
[xm(:)', fm]

func_obj(xm)



