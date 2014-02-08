function [x,fval,exitflag,output,population,score] = opti(nvars,lb,ub,PopulationSize_Data,Generations_Data)
% This is an auto generated M-file from Optimization Tool.

% Start with the default options
options = gaoptimset;
% Modify options setting
options = gaoptimset(options,'PopulationSize', PopulationSize_Data);
options = gaoptimset(options,'Generations', Generations_Data);
options = gaoptimset(options,'Display', 'off');
options = gaoptimset(options,'OutputFcns', { [] });
[x,fval,exitflag,output,population,score] = ...
ga(@hess_smith,nvars,[],[],[],[],lb,ub,[],options);
