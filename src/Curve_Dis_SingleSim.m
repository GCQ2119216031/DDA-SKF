%% ROC curve
load Gragh_Dis_singleSim
plot(SKF_FPR , SKF_TPR,'Color','#FF4500','LineWidth',1.5,'LineStyle','-');
hold on
plot(InterSim_FPR ,InterSim_TPR,'Color','#4169E1','LineWidth',1.5,'LineStyle','-');
hold on
plot(SemSim_FPR ,SemSim_TPR,'Color','#00C957','LineWidth',1.5,'LineStyle','-');
grid on
h2 = legend({'SKF AUROC = 0.828','Association Similarity AUROC = 0.707','Semantic Similarity AUROC = 0.741'});
set(h2,'FontName','Times New Roman','FontSize',10)
title('Disease Subspace ROC Curve','FontName','Times New Roman','FontSize',12)
xlabel('False Positive Rate','FontName','Times New Roman','FontSize',12)
ylabel('True Positive Rate','FontName','Times New Roman','FontSize',12)

%% PR curve
load Gragh_Dis_singleSim
plot(SKF_TPR ,SKF_Precision,'Color','#FF4500','LineWidth',1.5,'LineStyle','-');
hold on
plot(InterSim_TPR ,InterSim_Precision,'Color','#4169E1','LineWidth',1.5,'LineStyle','-');
hold on
plot(SemSim_TPR ,SemSim_Precision,'Color','#00C957','LineWidth',1.5,'LineStyle','-');
grid on
h2 = legend({'SKF AUPR = 0.359','Association Similarity AUPR = 0.214','Semantic Similarity AUPR = 0.108'});
set(h2,'FontName','Times New Roman','FontSize',10)
title('Disease Subspace PR Curve','FontName','Times New Roman','FontSize',12)
xlabel('Recall','FontName','Times New Roman','FontSize',12)
ylabel('Precision','FontName','Times New Roman','FontSize',12)