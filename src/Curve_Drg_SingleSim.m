%% ROC curve
load Gragh_Drug_singleSim
plot(SKF_FPR , SKF_TPR,'Color','#FF4500','LineWidth',1.5,'LineStyle','-');
hold on
plot(InterSim_FPR ,InterSim_TPR,'Color','#4169E1','LineWidth',1.5,'LineStyle','-');
hold on
plot(GoSim_FPR ,GoSim_TPR,'Color','#FFA500','LineWidth',1.5,'LineStyle','-');
hold on
plot(ChemSim_FPR ,ChemSim_TPR,'Color','#00C957','LineWidth',1.5,'LineStyle','-');
grid on
h2 = legend({'SKF AUROC = 0.888','Association Similarity AUROC = 0.812','Functional Similarity AUROC = 0.804','Chemical Similarity AUROC = 0.780'});
set(h2,'FontName','Times New Roman','FontSize',10)
title('Drug Subspace ROC Curve','FontName','Times New Roman','FontSize',12)
xlabel('False Positive Rate','FontName','Times New Roman','FontSize',12)
ylabel('True Positive Rate','FontName','Times New Roman','FontSize',12)

%% PR curve
load Gragh_Drug_singleSim
plot(SKF_TPR ,SKF_Precision,'Color','#FF4500','LineWidth',1.5,'LineStyle','-');
hold on
plot(InterSim_TPR ,InterSim_Precision,'Color','#4169E1','LineWidth',1.5,'LineStyle','-');
hold on
plot(GoSim_TPR ,GoSim_Precision,'Color','#FFA500','LineWidth',1.5,'LineStyle','-');
hold on
plot(ChemSim_TPR ,ChemSim_Precision,'Color','#00C957','LineWidth',1.5,'LineStyle','-');
grid on
h2 = legend({'SKF AUPR = 0.458','Association Similarity AUPR = 0.163','Functional Similarity AUPR = 0.111','Chemical Similarity AUPR = 0.092'});
set(h2,'FontName','Times New Roman','FontSize',10)
title('Drug Subspace PR Curve','FontName','Times New Roman','FontSize',12)
xlabel('Recall','FontName','Times New Roman','FontSize',12)
ylabel('Precision','FontName','Times New Roman','FontSize',12)