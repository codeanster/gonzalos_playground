# -*- coding: utf-8 -*-
"""Functions.ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/1ljUelj9eeTfolcrInNvwnTgIlZe5r-ES
"""

from IPython.display import clear_output
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import confusion_matrix, accuracy_score


def progress_bar(progress):
    """ Creates a progress bar of an execution """
    bar_length = 20
    if isinstance(progress, int):
        progress = float(progress)
    if not isinstance(progress, float):
        progress = 0
    if progress < 0:
        progress = 0
    if progress >= 1:
        progress = 1
        
    block = int(round(bar_length * progress))
    clear_output(wait = True)
    text = "Progress: [{0}] {1:.1f}%".format( "#" * block + "-" * (bar_length - block), progress * 100)
    print(text)


#for right now, assume data has been preprocessed
#in future, make pop up for user to choose, typing it out is barbaric 
def try_classification_models(X_train,y_train,X_test,y_test):
  """ Use all ML classifications on a problem and compare the accuracy 
      Inputs: X - test set (N of samples,p1,p2,pn+1)
              y - test labels (N of samples, 1)
              """
  assert X_train.shape[0] == y_train.shape[0], "Numer of samples has to be equal"
  assert (y_train.shape[-1] == 1) ,"y needs to be in shape of (N,1) explicitly"
  

  # Logistic_Regression 
  
  #make linear regression object
  classifier = LogisticRegression(random_state = 0)
  #fit
  classifier.fit(X_train, y_train)
  #predict
  classifier.predict(X_test)
  #accuracy measurements
  cm = confusion_matrix(y_test, y_pred)
  accuracy_score(y_test, y_pred)

  array = [cm,accuracy_score]
  return array

