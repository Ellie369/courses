#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 22 12:43:55 2020

@author: lily
"""
import matplotlib.pyplot as plt

history ={'loss': [0.6135852062225342,
  0.4178453644330303,
  0.36323275737166405,
  0.3319740912725528,
  0.3051732963959376,
  0.2853314482957125,
  0.26925592600057524,
  0.2600263031937182,
  0.24522654890343548,
  0.23199603473643463],
 'accuracy': [0.7775,
  0.8481333,
  0.8655,
  0.876,
  0.8862,
  0.89136666,
  0.8981,
  0.9025,
  0.9073667,
  0.9109667],
 'val_loss': [0.4725123842557271,
  0.4300117138773203,
  0.4162633110665613,
  0.40738544385466313,
  0.4128218510498603,
  0.4193886261847284,
  0.42967324300358695,
  0.4365849799166123,
  0.442405345539252,
  0.433763043437567],
 'val_accuracy': [0.841,
  0.853,
  0.851,
  0.86,
  0.862,
  0.857,
  0.867,
  0.861,
  0.857,
  0.867]}

plt.plot(history['accuracy'], label='train_accuracy')
plt.plot(history['val_accuracy'], label = 'test_accuracy')
plt.xlabel('Epoch')
plt.ylabel('Accuracy')
plt.ylim([0.5, 1])
plt.legend(loc='lower right')
plt.title('Training and Testing Accuracy of CNN')
plt.savefig('accuracy_cnn')