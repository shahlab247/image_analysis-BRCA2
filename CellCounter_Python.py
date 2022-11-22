import numpy as np
import time, os, sys
import matplotlib.pyplot as plt
import matplotlib as mpl
import fnmatch 
#%matplotlib inline
#mpl.rcParams['figure.dpi'] = 300
from cellpose import utils, io
import torch 
print(torch.cuda.get_device_name(0))
#%%
#importing tiff files

files = []

#enter the folder path which contains the tif files. 
input_dir = os.path.dirname()
for f in os.listdir(input_dir):
    filename = os.path.join(input_dir, f)
    files.append(filename)

#Selecting GFP files only
GFP_files=[]
for file in files:
    if fnmatch.fnmatch(file, '*GFP_RG*.tif'):
        GFP_files.append(file)

#to check if the right channel is imported. 
# view 1 image
img = io.imread(GFP_files[0])
plt.figure(figsize=(2,2))
plt.imshow(img)
plt.axis('off')
plt.show()

#%%
#RUN cellpose

from cellpose import models, io

# DEFINE CELLPOSE MODEL
# model_type='cyto' or model_type='nuclei'
#model = models.Cellpose(gpu=False, model_type='cyto')
model = models.CellposeModel(gpu=True, model_type='PCC_20X',device=torch.device("cuda:0"))
#model = models.CellposeModel(gpu=True, model_type='MTE firsttry')
# define CHANNELS to run segementation on
# grayscale=0, R=1, G=2, B=3`
# channels = [cytoplasm, nucleus]
# if NUCLEUS channel does not exist, set the second channel to 0
# channels = [0,0]
# IF ALL YOUR IMAGES ARE THE SAME TYPE, you can give a list with 2 elements
# channels = [0,0] # IF YOU HAVE GRAYSCALE
# channels = [2,3] # IF YOU HAVE G=cytoplasm and B=nucleus
# channels = [2,1] # IF YOU HAVE G=cytoplasm and R=nucleus

# or if you have different types of channels in each image
channels = [2,0]

#if diameter is set to None, the size of the cells is estimated on a per image basis
# you can set the average cell `diameter` in pixels yourself (recommended) 
# diameter can be a list or a single number for all images
                              
diam_labels = model.diam_labels.copy()                             

#defining a null array for saving filenames and number of cells. 
n_cell= [[2, 0]]*len(GFP_files) 

i=-1
# or in a loop
for filename in GFP_files:
    print(filename)
    img = io.imread(filename)
    masks, flows, diams = model.eval(img, diameter = diam_labels, channels=channels)
    i= i+ 1
    n_cell[i]= [filename,np.max(masks)]
    print('mask complete')
    