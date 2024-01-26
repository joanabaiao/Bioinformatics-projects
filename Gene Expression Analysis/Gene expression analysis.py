############ Gene Expression Analysis ##############

import numpy as np
from sklearn.model_selection import train_test_split
from sklearn import svm
from sklearn.neighbors import KNeighborsClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import confusion_matrix
from sklearn.metrics import roc_curve, auc
from sklearn.metrics import plot_confusion_matrix
import matplotlib.pyplot as plt

# IMPORT PROCESSED DATA
labels = np.loadtxt('E-GEOD-29044_Labels_proc.csv')
features = np.loadtxt('E-GEOD-29044_Features_proc.csv', delimiter = ';')

# DIVISION INTO TEST AND TRAINING SETS
X_train, X_test, y_train, y_test = train_test_split(features, labels, test_size = 0.3, random_state = 0)


############# SVM (Support vector machines) ############

clf_SVM =  svm.SVC(kernel = 'linear', C = 1).fit(X_train, y_train)
score_SVM = clf_SVM.score(X_test, y_test)

y_pred_SVM = clf_SVM.predict(X_test)

cm_SVM = confusion_matrix(y_test, y_pred_SVM)
TN_SVM, FP_SVM, FN_SVM, TP_SVM = cm_SVM.ravel()

precision_SVM = TP_SVM / (TP_SVM + FP_SVM)
sensitivity_SVM = TP_SVM / (TP_SVM + FN_SVM)
specificity_SVM = TN_SVM / (TN_SVM + FP_SVM)
f1_SVM = (2 * precision_SVM * sensitivity_SVM)/(precision_SVM + sensitivity_SVM)

y_score_SVM = clf_SVM.decision_function(X_test)
FP_rate_SVM, TP_rate_SVM, thresholds_SVM = roc_curve(y_test, y_score_SVM)
roc_auc_SVM = auc(FP_rate_SVM, TP_rate_SVM) 

print("---------- SVM CLASSIFIER ----------\n")
print("Confusion Matrix:\n", cm_SVM)
print('\n* Accuracy:', round(score_SVM*100, 3), '%')
print('* Precision:', round(precision_SVM*100, 3), '%')
print('* Sensitivity:', round(sensitivity_SVM*100, 3), '%')
print('* Specificity:', round(specificity_SVM*100, 3), '%')
print("* F1 Score:", round(f1_SVM*100, 3), '%')
print("* AUC:", round(roc_auc_SVM, 3), '\n')

plot_confusion_matrix(clf_SVM, X_test, y_test) 
plt.title("Confusion Matrix - SVM")
plt.show()

plt.plot(FP_rate_SVM, TP_rate_SVM, 'b', label = "AUC = %0.3f" % roc_auc_SVM)
plt.plot([0, 1], [0, 1], 'r--')
plt.legend(loc = "lower right")
plt.xlim([-0.1, 1.1])
plt.ylim([-0.1, 1.1])
plt.title("ROC Curve - SVM")
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.show()


############### KNN (K-nearest neighbors) #############

# Find the ideal K value:
error_rate = []
for i in range(1,30):
    knn = KNeighborsClassifier(n_neighbors = i).fit(X_train,y_train)
    pred = knn.predict(X_test)
    error_rate.append(np.mean(pred != y_test))

plt.plot(range(1, 30), error_rate, linestyle='dashed', marker='o')
plt.title('Error Rate vs. Value of k')
plt.xlabel('k')
plt.ylabel('Error Rate')
k = error_rate.index(min(error_rate)) + 1

# Classification + Performance Metrics
clf_KNN = KNeighborsClassifier(n_neighbors = k).fit(X_train, y_train)
score_KNN = clf_KNN.score(X_test, y_test)

y_pred_KNN = clf_KNN.predict(X_test)

cm_KNN = confusion_matrix(y_test, y_pred_KNN)
TN_KNN, FP_KNN, FN_KNN, TP_KNN = cm_KNN.ravel()
precision_KNN = TP_KNN / (TP_KNN + FP_KNN)
sensitivity_KNN = TP_KNN / (TP_KNN + FN_KNN)
specificity_KNN = TN_KNN / (TN_KNN + FP_KNN)
f1_KNN = (2 * precision_KNN * sensitivity_KNN)/(precision_KNN + sensitivity_KNN)

y_score_KNN = clf_KNN.predict_proba(X_test) #duas colunas (classe negativa e outra para a classe positiva)
FP_rate_KNN, TP_rate_KNN, thresholds_KNN = roc_curve(y_test, y_score_KNN[:,1]) #toda a coluna da classe positiva
roc_auc_KNN = auc(FP_rate_KNN, TP_rate_KNN)


print("\n---------- KNN CLASSIFIER ----------\n")
print("Minimum Error:", round(min(error_rate), 3), "for K =", k)
print("\nConfusion Matrix:\n", cm_KNN)
print("\n* Accuracy:", round(score_KNN*100, 3), '%')
print('* Precision:', round(precision_KNN*100, 3), '%')
print('* Sensitivity:', round(sensitivity_KNN*100, 3), '%')
print('* Specificity:', round(specificity_KNN*100, 3), '%')
print("* F1 Score:", round(f1_KNN*100, 3), '%')
print("* AUC:", round(roc_auc_KNN, 3), '\n')


plot_confusion_matrix(clf_KNN, X_test, y_test)  
plt.title("Confusion Matrix - KNN")
plt.show()

plt.plot(FP_rate_KNN, TP_rate_KNN, 'b', label="AUC = %0.3f" % roc_auc_KNN)
plt.plot([0, 1], [0, 1], 'r--')
plt.legend(loc="lower right")
plt.xlim([-0.1, 1.1])
plt.ylim([-0.1, 1.1])
plt.title("ROC Curve - KNN")
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.show()


##################### CLASSIFICADOR RANDOM FOREST ###################

clf_RF = RandomForestClassifier(n_estimators = 20, random_state = 0).fit(X_train, y_train)
score_RF = clf_RF.score(X_test, y_test)

y_pred_RF = clf_RF.predict(X_test)

cm_RF = confusion_matrix(y_test, y_pred_RF)
TN_RF, FP_RF, FN_RF, TP_RF = cm_RF.ravel()
precision_RF = TP_RF / (TP_RF + FP_RF)
sensitivity_RF = TP_RF / (TP_RF + FN_RF)
specificity_RF = TN_RF / (TN_RF + FP_RF)
f1_RF = (2 * precision_RF * sensitivity_RF)/(precision_RF + sensitivity_RF)

y_score_RF = clf_RF.predict_proba(X_test) 
FP_rate_RF, TP_rate_RF, thresholds_RF = roc_curve(y_test, y_score_RF[:,1]) 
roc_auc_RF = auc(FP_rate_RF, TP_rate_RF)

print("\n---------- RANDOM FOREST CLASSIFIER ----------\n")
print("Confusion Matrix:\n", cm_RF)
print('\n* Accuracy: ', round(score_RF*100, 3), '%')
print('* Precision:', round(precision_RF*100, 3), '%')
print('* Sensitivity:', round(sensitivity_RF*100, 3), '%')
print('* Specificity:', round(specificity_RF*100, 3), '%')
print("* F1 Score:", round(f1_RF*100, 3), '%')
print("* AUC:", round(roc_auc_RF, 3), '\n')

plot_confusion_matrix(clf_RF, X_test, y_test)  
plt.title("Confusion Matrix - RF")
plt.show()

plt.plot(FP_rate_RF, TP_rate_RF, 'b', label="AUC = %0.3f" % roc_auc_RF)
plt.plot([0, 1], [0, 1], 'r--')
plt.legend(loc="lower right")
plt.xlim([-0.1, 1.1])
plt.ylim([-0.1, 1.1])
plt.title("ROC Curve - RF")
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.show()