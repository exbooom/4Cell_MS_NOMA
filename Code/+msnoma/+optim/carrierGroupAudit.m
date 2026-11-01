function audit = carrierGroupAudit(num)
%CARRIERGROUPAUDIT Expose Eq24 overrun; never use truncated groups to solve.
width=2*num.G-1;
audit.printedWidth=width;
audit.printedFirst=(0:num.M-1)'*width+1;
audit.printedLast=(1:num.M)'*width;
audit.printedCarrierCount=num.M*width;
audit.actualCarrierCount=num.N;
audit.outOfBandCount=max(0,num.M*width-num.N);
audit.uncoveredCount=max(0,num.N-num.M*width);
audit.isExactPartition=num.M*width==num.N;
audit.solverUsesPrintedGroups=false;
audit.candidatePolicy='full_A_columns_all_communication_receivers_all_incoming_hearability_rows';
end
