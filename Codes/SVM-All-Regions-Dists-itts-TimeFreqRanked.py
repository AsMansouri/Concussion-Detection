#Import scikit-learn dataset library
from sklearn.model_selection import train_test_split
from sklearn.metrics import confusion_matrix
from sklearn import svm
from sklearn import metrics
from scipy.io import loadmat
#import tensorflow as tf
import numpy as np
import xlsxwriter 
import glob
import os.path
import os

Regions = {'R-Frontal':np.array([1,2,3,4,5,6,10,11,12,13,14,15,18,19,20,21,25,26,222,223,224])-1,
           'L-Frontal':np.array([21,22,23,26,27,28,29,32,33,34,35,36,37,38,39,40,46,47,48,54])-1,
           'Frontal':np.array([1,2,3,4,5,6,10,11,12,13,14,15,18,19,20,21,25,26,222,223,224,22,23,26,27,28,29,32,33,34,35,36,37,38,39,40,46,47,48,54])-1,
           'R-Temporal':np.array([170,171,172,178,179,180,181,190,191,192,193,194,202,203,204,210,211,212,220,221])-1,
           'L-Temporal':np.array([55,56,57,61,62,63,64,68,69,70,71,74,75,76,83,84,85,94,95,96])-1,
           'Temporal':np.array([170,171,172,178,179,180,181,190,191,192,193,194,202,203,204,210,211,212,220,221,55,56,57,61,62,63,64,68,69,70,71,74,75,76,83,84,85,94,95,96])-1,
           'R-Central':np.array([7,8,81,90,130,131,132,142,143,144,154,155,163,164,173,182,183,184,185,186,195,196,197,198,205,206,207,213,214,215])-1,
           'L-Central':np.array([8,9,16,17,24,30,41,42,43,44,45,49,50,51,52,53,58,59,60,65,66,72,77,78,79,80,81,88,89,90])-1,
           'Central':np.array([7,8,81,90,130,131,132,142,143,144,154,155,163,164,173,182,183,184,185,186,195,196,197,198,205,206,207,213,214,215,9,16,17,24,30,41,42,43,44,45,49,50,51,52,53,58,59,60,65,66,72,77,78,79,80,88,89])-1,
           'Parietal':np.array([85,86,87,97,98,99,100,101,108,109,110,118,119,127,128,129,140,141,151,152,153,161,162,171])-1,
           'Occipital':np.array([106,107,108,114,115,116,117,123,124,125,126,136,137,138,139,147,148,149,150,151,158,159,160,168,169])-1,
           'Right':np.array([1,2,3,4,5,6,7,8,10,11,12,13,14,15,18,19,20,21,25,26,81,90,101,119,126,127,128,129,130,131,132,137,138,139,140,141,142,143,144,148,149,150,151,152,153,154,155,158,159,160,161,162,163,164,168,169,170,171,172,173,178,179,180,181,182,183,184,185,186,190,191,192,193,194,195,196,197,198,202,203,204,205,206,207,210,211,212,213,214,215,220,221,222,223,224])-1,
           'Left':np.array([8,9,16,17,21,22,23,24,26,27,28,29,30,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,68,69,70,71,72,74,75,76,77,78,79,80,81,83,84,85,86,87,88,89,90,94,95,96,97,98,99,100,101,106,107,108,109,110,114,115,116,117,118,123,124,125,126,136,137])-1,
           'All':np.array([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,68,69,70,71,72,74,75,76,77,78,79,80,81,83,84,85,86,87,88,89,90,94,95,96,97,98,99,100,101,105,106,107,108,109,110,114,115,116,117,118,119,123,124,125,126,127,128,129,130,131,132,136,137,138,139,140,141,142,143,144,147,148,149,150,151,152,153,154,155,158,159,160,161,162,163,164,168,169,170,171,172,173,177, 178,179,180,181,182,183,184,185,186,190,191,192,193,194,195,196,197,198,202,203,204,205,206,207,210,211,212,213,214,215,220,221,222,223,224])-1,
           'All 256 ch.':range(1,256)}

print("Reading Dists Datasets...")
files = glob.glob("../Datasets_TimeFreqRanked/*.mat")
if not(os.path.exists('../Results/')):
    os.mkdir('../Results/')
if not(os.path.exists('../Results/ResultsTimeFreqRanked/')):
    os.mkdir('../Results/ResultsTimeFreqRanked/')

for filepath in files:
    #print(filepath)
    tmp = filepath.split("PrePost_",1)[1][:-4]
    OutFile = ('../Results/ResultsTimeFreqRanked/SVM_Dist_' + tmp + '.txt')
    #workbook = xlsxwriter.Workbook('Results/Dist_' + tmp + '.xlsx') 
    #worksheet = workbook.add_worksheet(tmp)

    if not(os.path.exists(OutFile)):

        f= open(OutFile,"w")
        f.close()

        with open(OutFile, 'a') as results_data:
                results_data.write('\n{}\t{}\t{}\t{}'.format('File Name', '#itt', 'Regions', 'ACC'))
        """
        with open(OutFile, 'a') as results_data:
                results_data.write('\n{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format('File Name', '#itt', 'Regions', 'ACC', 'Precision', 'TPR', 'TP', 'TN', 'FN', 'FP'))
        """
        
        x_test = loadmat(filepath)['Dataset_Test']
        x_train = loadmat(filepath)['Dataset_Train']
        y_test = loadmat(filepath)['Labels_Test']
        y_train = loadmat(filepath)['Labels_Train']

        y_test = y_test.reshape(len(y_test))
        y_train = y_train.reshape(len(y_train))

        #print(y_test.shape, y_train.shape)
        #print(x_test.shape, x_train.shape)

        for itt in range(1):
            try:
                print(filepath + " - itt = " + str(itt))

                for i in range(len(Regions.keys())):
                    if x_train.shape[1] > 256:
                        tmpregionsCH = np.append(np.array(list(Regions.values())[i]),np.array(list(Regions.values())[i])+256)
                        x_train_tmp_tmp = x_train[:,tmpregionsCH,:]
                        x_train_tmp = x_train_tmp_tmp[:,:,list(Regions.values())[i]]
                        x_test_tmp_tmp = x_test[:,tmpregionsCH,:]
                        x_test_tmp = x_test_tmp_tmp[:,:,list(Regions.values())[i]]
                    else:
                        x_train_tmp_tmp = x_train[:,list(Regions.values())[i],:]
                        x_train_tmp = x_train_tmp_tmp[:,:,list(Regions.values())[i]]
                        x_test_tmp_tmp = x_test[:,list(Regions.values())[i],:]
                        x_test_tmp = x_test_tmp_tmp[:,:,list(Regions.values())[i]]
                
                    x_train1 = x_train_tmp.reshape(x_train_tmp.shape[0],x_train_tmp.shape[1]*x_train_tmp.shape[2])        
                    x_test1 = x_test_tmp.reshape(x_test_tmp.shape[0],x_test_tmp.shape[1]*x_test_tmp.shape[2])
                    x_train1[np.where(np.isnan(x_train1))] = 0 
                    x_test1[np.where(np.isnan(x_test1))] = 0 

                    #Create a svm Classifier
                    clf = svm.SVC(kernel='linear', decision_function_shape='ovr') # Linear Kernel
                    #Train the model using the training sets
                    clf.fit(x_train1, y_train)

                    y_pred = clf.predict(x_test1)
                    CF = confusion_matrix(y_test,y_pred)


                    with open(OutFile, 'a') as results_data:
                        results_data.write('\n{}\t{}\t{}\t{}'.format(filepath, str(itt), list(Regions.keys())[i], metrics.accuracy_score(y_test, y_pred)))
                        
                    """
                    with open(OutFile, 'a') as results_data:
                        results_data.write('\n{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(filepath, str(itt), list(Regions.keys())[i], 
                        metrics.accuracy_score(y_test, y_pred), metrics.precision_score(y_test, y_pred), metrics.recall_score(y_test, y_pred),
                        CF[0][0], CF[1][1], CF[0][1]), CF[1][0])
                    """

            except ValueError:
                print('Value Error!!')
            else:
                #print("Couldn't perform on" + print(filepath))
                #worksheet = workbook.add_worksheet(filepath[-34:-4])
                #worksheet.write(0, 0, 'Could not perform on this file')
                continue
    else:
        print(OutFile + " It exist...") 
        continue               
        