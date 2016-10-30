figure()
errorbar(MassFlow_TEV,1,MassFlowError_TEV,'-xg')
hold on
errorbar(MassFlow_CT,2,MassFlowError_CT,'-*r')
errorbar(MassFlow_HEV1,3,MassFlowError_HEV1,'-^b')
errorbar(MassFlow_HEV2,4,MassFlowError_HEV2,'-sm')
errorbar(MassFlow_HEV3,5,MassFlowError_HEV3,'-dc')
legend('TEV', 'CP', 'HEV1','HEV2','HEV3','Location','Southeast')
xlabel('Valve')
ylabel('Mass Flow Rate')
title('Mass Flow Rate vs Valve')
hold off

figure()
errorbar(W_comp_TEV,1,u_comp_work_TEV,'-xg')
hold on
errorbar(W_comp_CT,2,u_comp_work_CT,'-*r')
errorbar(W_comp_HEV1,3,u_comp_work_HEV1,'-^b')
errorbar(W_comp_HEV2,4,u_comp_work_HEV2,'-sm')
errorbar(W_comp_HEV3,5,u_comp_work_HEV3,'-dc')
legend('TEV', 'CP', 'HEV1','HEV2','HEV3','Location','Southeast')
xlabel('Valve')
ylabel('Compressor Work')
title('Compressor Work vs Valve')
hold off

figure()
errorbar(W_cooling_pow_TEV,1,u_cooling_pow_TEV,'-xg')
hold on
errorbar(W_cooling_pow_CT,2,u_cooling_pow_CT,'-*r')
errorbar(W_cooling_pow_HEV1,3,u_cooling_pow_HEV1,'-^b')
errorbar(W_cooling_pow_HEV2,4,u_cooling_pow_HEV2,'-sm')
errorbar(W_cooling_pow_HEV3,5,u_cooling_pow_HEV3,'-dc')
legend('TEV', 'CP', 'HEV1','HEV2','HEV3','Location','Southeast')
xlabel('Valve')
ylabel('Cooling Power')
title('Cooling Power vs Valve')
hold off


figure()
errorbar(COP_TEV,1,uCOP_TEV,'-xg')
hold on
errorbar(COP_CT,2,COP_CT,'-*r')
errorbar(COP_HEV1,3,COP_HEV1,'-^b')
errorbar(COP_HEV2,4,COP_HEV2,'-sm')
errorbar(COP_HEV3,5,COP_HEV3,'-dc')
legend('TEV', 'CP', 'HEV1','HEV2','HEV3','Location','Southeast')
xlabel('Valve')
ylabel('CoP')
title('CoP vs Valve')
hold off

