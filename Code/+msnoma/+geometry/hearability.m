function out = hearability(powerW,ch,spec)
%HEARABILITY All-pair received power rule; no absolute acquisition threshold.
received = ch.gP.*powerW;
[K,M] = size(received);
out.mask = false(K,M);
out.ratioToStrongestOther = zeros(K,M);
for k = 1:K
    other = max(received([1:k-1,k+1:K],:),[],1);
    active = received(k,:)>0;
    ratio = inf(1,M);
    positive = other>0;
    ratio(positive) = received(k,positive)./other(positive);
    ratio(~active) = 0;
    out.ratioToStrongestOther(k,:) = ratio;
    out.mask(k,:) = active & ratio>=spec.ratio*(1-spec.comparisonTolerance);
end
out.receivedPowerW = received;
out.count = sum(out.mask,1)';
out.enoughBs = out.count>=spec.minimumBs;
out.model = 'relative_power_rule_only_not_measured_detection_probability';
end
