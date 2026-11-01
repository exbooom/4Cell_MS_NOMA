function exportSubgradient(run,runDir)
%EXPORTSUBGRADIENT Preserve every inner call, including rejected outer trials.
index=find(strcmp({run.algorithms.id},'cp4a_subgradient_candidate'),1);
if isempty(index), return; end
s=run.algorithms(index).solution;
if ~isfield(s,'innerTrace'), return; end
count=numel(s.innerTrace);
data=zeros(count,8+run.num.K); reason=cell(count,1);
allHistory=zeros(0,9);
for i=1:count
    item=s.innerTrace{i}; st=item.state;
    data(i,:)=[i,st.iterations,st.evaluations,double(st.converged), ...
        st.maxVariableViolation,st.innerComplementarity,st.innerRelativeComplementarity, ...
        -st.negativeDual,item.budgetPrice'];
    reason{i}=st.reason;
    allHistory=[allHistory;[repmat(i,size(st.history,1),1),st.history]]; %#ok<AGROW>
end
names=[{'innerCall','iterations','functionEvaluations','converged', ...
    'maxVariableViolation','maxComplementarity','relativeAbsoluteComplementarity', ...
    'scaledDualBound'},cellstr(compose('budgetPriceBS%d',1:run.num.K))];
t=array2table(data,'VariableNames',names); t.reason=reason;
writetable(t,fullfile(runDir,'subgradient_inner_calls.csv'));
writetable(array2table(allHistory,'VariableNames',[{'innerCall'}, ...
    s.innerTrace{1}.state.historyColumns]),fullfile(runDir,'subgradient_inner_history.csv'));
prices=array2table([s.history(:,1),s.outerPriceHistory], ...
    'VariableNames',[{'outerIteration'},cellstr(compose('budgetPriceBS%d',1:run.num.K))]);
writetable(prices,fullfile(runDir,'subgradient_outer_prices.csv'));
end
