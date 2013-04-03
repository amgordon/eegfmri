function r = sliding_window(nx, nwind, noverlap)

r = bsxfun(@plus, (1:nwind)', 1+(0:(fix((nx-noverlap)/(nwind-noverlap))-1))*(nwind-noverlap))-1;