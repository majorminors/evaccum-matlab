function [bestpar,bestval,BIC]=fitparams_refine_template(errfunc,mod_feature,data2fit,randiter,nosession,startpar,parange)

nopara=length(parange);     % number of parameters
bestpar = zeros (nosession, nopara);    % best fit parameters
bestval = Inf (nosession, 1) ;          % best error functions
BIC=Inf(nosession,1); %start big for optimisation
goalstat = cell(1, length(data2fit));    % data to fit

if length(data2fit)==1
    if isstruct(data2fit)
        goalstat=data2fit;
    else
        goalstat=data_stats(data2fit);
    end
else
    if isstruct(data2fit{1})
        goalstat=data2fit;
    else
        for j=1:length(data2fit)
            goalstat{j} = data_stats(data2fit{j});
        end
    end
end


% do the optimization
for iter = 1:nosession
    disp ('------------------------ ');
    disp(sprintf ('Starting optimization session number: %d', iter));
    if randiter == 0
        bestpar (iter,:) = startpar;
    else
        disp ('Searching for starting point of optimization');
        bestval(iter)=nan;
        while isnan(bestval(iter))
        for i = 1:randiter
            if isempty(startpar)
                param = ((parange(:,2)-parange(:,1)) .* rand (1, nopara)' + parange(:,1))';
                param(2)=rand*param(1);
            else
                minrange = max (parange(:,1)'./startpar, 0.0001);
                maxrange = min (parange(:,2)'./startpar, 2);
                param = (maxrange-minrange) .* rand (1, nopara) + minrange;
                param(2)=rand(1,1)*param(1); % C0 is smaller than threshold
            end
            erp = feval(errfunc,param, mod_feature, goalstat);
            if erp < bestval(iter) || i == 1
                bestpar(iter,:) = param;
                bestval(iter) = erp;
            end
        end
        end
    end
    disp (['Optimizing parameters, the ' num2str(iter) ' iterations, starting fitness: ' num2str(bestval(iter))]);
    bestpar(iter,:)=fitparams_template(errfunc,bestpar(iter,:),mod_feature,goalstat,parange);
    [bestval(iter),BIC(iter)]=feval(errfunc,bestpar(iter,:), mod_feature, goalstat);
%     clear foo2 foo3 foo4 foo5 foo6 foo7;
    disp (['Optimizing parameters finished, final fitness: ' num2str(bestval(iter))]);
end