function out = varianceThreeTerm(powerW,ch,num,factor)
%VARIANCETHREETERM Eq9: evaluate omitted positioning interference post-solve.
out = msnoma.ranging.varianceTwoTerm(powerW,factor);
received = ch.gP.*powerW;
otherReceivedW = max(0,sum(received,1)-received);
extra = factor.aPerS*num.TpS/(2*num.BfeHz^2)*otherReceivedW./received;
extra(received==0) = Inf;
out.positioningDelayVarianceS2 = extra;
out.delayVarianceS2 = out.delayVarianceS2+extra;
out.rangeVarianceM2 = out.delayVarianceS2*num.c0Mps^2;
out.rangeStdM = sqrt(out.rangeVarianceM2);
out.model = 'paper_Eq9_three_term_post_evaluation_not_fixed_R_optimization';
end
