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
                   train_x[num_image] = cv2.resize(temp.pixel_array, (std_width, std_height), interpolation = cv2.INTER_AREA)
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
    train_x_sharpened[i] = sharpen_image(train_x[i])

# In[6]:
# find index of black image, whose first row's sum is zero
def find_black_index(image):
    _list = []
    for i in range(0, len(image)):
        if (sum(image[i, 0]) > ):
            _list.append(i)
    return _list

# remove index of lists
def remove_index(image, _list):
    image_index = [i for i in range(len(image))]
    image_index = set(image_index)
    black_index = set(_list)
    index = list(image_index - black_index)
    new_image = [image[i] for i in index]
    return new_image

black_index = find_black_index(train_x)
train_x_re = remove_index(train_x, black_index)
train_y_re = remove_index(train_y, black_index)


# In[7]:
def get_FileSize(filePath):
    # filePath = unicode(filePath, 'utf8')
    fsize = os.path.getsize(filePath)
    fsize = fsize / float(1024 * 1024)
    return round(fsize, 2)


# In[8]:
train_y_cala = pd.read_csv(***YOUR PATH***)
train_y_cala = [train_y_cala.assessment, train_y_cala.pathology]

std_height = 299
std_width = 299
train_x_cala = np.zeros(shape = (10, std_height, std_width))
num_image_cala = 0
train_x_cala_ori_size = np.zeros(shape = (10, 2))
train_x_index = np.zeros(shape = (2 * 10))
train_x_missing = []
train_x_total_num = 0
# rootdit: CBIS-DDSM is default  or other folder where you save images
def list_all_files_by_size(rootdir):
    global train_x_cala, num_image_cala, train_x_missing, train_x_total_num, train_x_index
    _files = []
    list = os.listdir(rootdir) # list all dir and file
    for i in range(0, len(list)):
        path = os.path.join(rootdir,list[i])
        if os.path.isdir(path):
            _files.extend(list_all_files_by_size(path))
        if os.path.isfile(path):
            # fsize = os.path.getsize(path)
            
            # read dcm image
            temp = pydicom.dcmread(path)
            if sum(sum(temp.pixel_array)) > 1000000:
                # resize and add new pixel array to train_x
                # interpolation is parameter of different methods to resize image
                train_x_cala_ori_size[num_image_cala] = [len(temp.pixel_array), len(temp.pixel_array[0])]
                train_x_cala[num_image_cala] = cv2.resize(temp.pixel_array, (std_width, std_height), interpolation = cv2.INTER_AREA)
                num_image_cala += 1
                train_x_index[train_x_total_num] = 1
                # continue
            
            _files.append(path)
            train_x_total_num += 1

    return _files

def find_missing(train_x_index, train_x_missing):
    for i in range(0, len(train_x_index)):
        if (i > 0) & (i % 2 == 1) & (train_x_index[i - 1] == 0) & (train_x_index[i] == 0):
            train_x_missing.append(int(i / 2))

def remove_index(image, _list):
    image_index = [i for i in range(len(image))]
    image_index = set(image_index)
    black_index = set(_list)
    index = list(image_index - black_index)
    new_image = [image[i] for i in index]
    return new_image

#_list = list_all_files_by_size(***YOUR PATH***)
_list = list_all_files_by_size('D:\\PR Project\\Project\\data\\CALA-Train')

find_missing(train_x_index, train_x_missing)

train_y_cala = remove_index(train_y_cala, train_x_missing)

# type your image folder
# _file_list = list_all_files('D:\\PR Project\\Project\\data\\CBIS-DDSM')
# after this, train_x should be already preprocessed into 3-D array in same shape


# %%
