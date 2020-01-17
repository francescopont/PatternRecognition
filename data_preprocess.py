# In[1]:
import pandas as pd
import numpy as np
import os
import pydicom
import cv2

# In[2]:
# https://wiki.cancerimagingarchive.net/display/Public/CBIS-DDSM#5e40bd1f79d64f04b40cac57ceca9272
# use corresponding image and csv file
train_data = pd.read_csv('D:\\PR Project\\Project\\data\\mass_case_description_train_set.csv')
train_y = [train_data.assessment, train_data.pathology]

# In[3]:
# standrad size
# change h and w to resize image if reach memory limitation
std_height = 300
std_width = 300

# initial train_x array of size: image number * std_h * std_w
train_x = np.zeros(shape = (len(train_y[0]), std_height, std_width))
num_image = 0
train_x_ori_size = np.zeros(shape = (len(train_y[0]), 2))

# rootdit: CBIS-DDSM is default  or other folder where you save images
def list_all_files(rootdir, num_image):
    global train_x, num_image
    _files = []
    list = os.listdir(rootdir) # list all dir and file
    for i in range(0, len(list)):
         path = os.path.join(rootdir,list[i])
         if os.path.isdir(path):
              _files.extend(list_all_files(path))
         if os.path.isfile(path):
              # read dcm image
              temp = pydicom.dcmread(path)

              # use ONLY for distinguish cropped image and ROI image
              # because they're in same floder
              if (len(temp.pixel_array) < 1000):
                   # resize and add new pixel array to train_x
                   # interpolation is parameter of different methods to resize image
                   train_x_ori_size[num_image] = [len(ds.pixel_array), len(ds.pixel_array[0])]
                   train_x[num_image] = cv2.resize(temp.pixel_array, (std_width, std_height), interpolation = cv2.INTER_CUBIC)
                   num_image += 1
                   # continue

              # resize and add new pixel array to train_x
              # interpolation is parameter of different methods to resize image
              # train_x[num_image] = cv2.resize(temp.pixel_array / norm_pixel, (std_width, std_height), interpolation = cv2.INTER_CUBIC)
              # num_image += 1
               
              _files.append(path)
    return _files
# type your image folder
_file_list = list_all_files('D:\\PR Project\\Project\\data\\CBIS-DDSM')
# after this, train_x should be already preprocessed into 3-D array in same shape

# In[4]:
# Normalization train_x
def normalization(data):
    _range = np.max(data) - np.min(data)
    return (data - np.min(data)) / _range
train_x = normalization(train_x)

# Binary train_y label
# to be added

# train_x[i] is pixel array for No. i image, whose shape is std_h * std_w, for now, 480 * 300
# train_y[0] is assessment result, and train_y[1] is pathology result
# for training No. i image, use train_x[i] as 2-D input, and use train_y[0][i] OR train_y[1][i] as output
# In[5]:
# sharpen image function using Laplacian operator
def sharpen_image(image):
    kernel = np.array([[0, -1, 0], [-1, 5, -1], [0, -1, 0]], np.float32) # kernal size
    dst = cv2.filter2D(image, -1, kernel = kernel)
    return dst

# sharpen image
train_x_sharpened = train_x
for i in range(0, len(train_x)):
    train_x[i] = sharpen_image(train_x[i])



# %%
