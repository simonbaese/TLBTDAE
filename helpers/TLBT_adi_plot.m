function TLBT_adi_plot(outB,outC)
%% Delivers plots of the normalized residual norms of the ADI iteration
%
% Inputs:
% outB, outC   output structs of LDLT ADI iteration (mess_lradi), must 
%              contain outB.res and outC.res
%
% Author: Simon Bäse (2021)
% Partly forked from DAE2 demo in M-M.E.S.S. toolbox by
% Saak, J. and Köhler, M. and Benner, P. and others
%
% See LICENSE.md

figure;
semilogy(outB.res,'LineWidth',1);
title('0 = BB^T + AXM^T + MXA^T');
set(gca,'XLimSpec','Tight');
xlabel('number of iterations');
ylabel('normalized residual norm');

figure;
semilogy(outC.res,'LineWidth',1);
title('0 = C^TC + A^TXE + E^TXA');
set(gca,'XLimSpec','Tight');
xlabel('number of iterations');
ylabel('normalized residual norm');    

end

