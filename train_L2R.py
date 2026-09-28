# This script imports training data and trains L2R MAE model
#
# Last edits: 09/28/2026 BG 

import matplotlib.pyplot as plt
import numpy as np
#np.long = int 
import pandas as pd
import tensorflow as tf
import os
import glob
import mat73
import gc
from scipy.io import loadmat
from sklearn.metrics import accuracy_score, precision_score, recall_score
from sklearn.model_selection import train_test_split
from tensorflow.keras import layers, losses
#from tensorflow.keras.datasets import fashion_mnist (fake data)
from tensorflow.keras.models import Model
from PIL import Image
import keras
from matplotlib.animation import FuncAnimation


# Paths
modelsavepath=r'/Volumes/kanarde/MasonBEAST/data/trained_models'
#'/Volumes/Elements/StormCHAZerz Data/Dec2023Noreaster_processed/1702827001820/MAE_data/trained_models/'
savepath=r'/Volumes/kanarde/MasonBEAST/data/testingtime1'
#'/Volumes/Elements/StormCHAZerz Data/Dec2023Noreaster_processed/1702827001820/MAE_data/testingtime1/'
transect_1702827001820_path=r'/Volumes/kanarde/MasonBEAST/data/DEMs/1702827001820/Transects/alongshore_transects.mat'

figpath1820=os.path.join(savepath,'1702827001820/Figures/')
savepath1820=os.path.join(savepath,'1702827001820/')

# Functions 
