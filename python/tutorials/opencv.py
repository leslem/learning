# https://towardsdatascience.com/computer-vision-for-beginners-part-1-7cca775f58ef

import cv2
import numpy as np
import matplotlib.pyplot as plt
import os


images_path = 'computer-vision-images'

## Reading in images and figuring out their color space

# The openCV default image color format is BGR?!?! so this image looks weird (orange sky)
img = cv2.imread(os.path.join(images_path, 'burano.jpg'))
plt.imshow(img)
plt.show()


# Convert it to RGB like a sane person.
img_rgb = cv2.cvtColor(img, cv2.COLOR_BGR2RGB)
plt.imshow(img_rgb)
plt.show()

# Conver to grayscale because why not
img_gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
plt.imshow(img_gray, cmap = 'gray')
plt.show()

# Plot the three channels of the image
fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(20, 20))

for i in range(0, 3):
    ax = axs[i]
    ax.imshow(img_rgb[:, :, i], cmap='gray')

plt.show()

# The image object is a numpy array of RGB pixel values.
img_rgb
dir(img_rgb)
img_rgb.shape
# 168 x 300 x 3
# This image is 168 pixels high x 300 pixels wide, and each pixel has RGB values


## Drawing on an image in matplotlib
img = cv2.cvtColor(cv2.imread(os.path.join(images_path, 'wall.jpg')), cv2.COLOR_BGR2RGB)
plt.imshow(img)
plt.show()

img_copy = img.copy()
cv2.rectangle(img_copy, pt1=(800, 470), pt2=(980, 530),
              color=(255, 0, 0), thickness=5)
plt.imshow(img_copy)
plt.show()

# https://towardsdatascience.com/computer-vision-for-beginners-part-2-29b3f9151874

## Filtering via CNN
img = cv2.cvtColor(cv2.imread(os.path.join(images_path, 'text.jpg')), cv2.COLOR_BGR2RGB)

# Plot image with different kernel sizes for blurring
kernel_sizes = [5, 11, 17]

fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(20, 20))
for ind, s in enumerate(kernel_sizes):
    img_blurred = cv2.blur(img, ksize=(s, s))  # Average blurring
    ax = axs[ind]
    ax.imshow(img_blurred)
    ax.axis('off')

plt.show()

# Bilateral filtering does better at blurring noise while also keeping edges sharp
fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(20, 20))
for ind, s in enumerate(kernel_sizes):
    img_blurred = cv2.bilateralFilter(img, s, sigmaSpace=75, sigmaColor=75)  # Bilateral filter
    ax = axs[ind]
    ax.imshow(img_blurred)
    ax.axis('off')

plt.show()

# Different ways of thresholding
_, thresh_0 = cv2.threshold(img, 127, 255, cv2.THRESH_BINARY)
_, thresh_1 = cv2.threshold(img, 127, 255, cv2.THRESH_BINARY_INV)
_, thresh_2 = cv2.threshold(img, 127, 255, cv2.THRESH_TOZERO)
_, thresh_3 = cv2.threshold(img, 127, 255, cv2.THRESH_TOZERO_INV)
_, thresh_4 = cv2.threshold(img, 127, 255, cv2.THRESH_TRUNC)

# Plot the images
images = [img, thresh_0, thresh_1, thresh_2, thresh_3, thresh_4]
fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(13, 13))
for ind, p in enumerate(images):
    ax = axs[ind // 3, ind % 3]
    ax.imshow(p)

# Adaptive Thresholding
_, thresh_binary = cv2.threshold(img, thresh=127, maxval=255, type=cv2.THRESH_BINARY)
adap_mean_2 = cv2.adaptiveThreshold(img, 255, 
                                    cv2.ADAPTIVE_THRESH_MEAN_C, 
                                    cv2.THRESH_BINARY, 7, 2)
adap_mean_2_inv = cv2.adaptiveThreshold(img, 255, 
                                        cv2.ADAPTIVE_THRESH_MEAN_C, 
                                        cv2.THRESH_BINARY_INV, 7, 2)
adap_mean_8 = cv2.adaptiveThreshold(img, 255, 
                                    cv2.ADAPTIVE_THRESH_MEAN_C, 
                                    cv2.THRESH_BINARY, 7, 8)
adap_gaussian_8 = cv2.adaptiveThreshold(img, 255, 
                                        cv2.ADAPTIVE_THRESH_GAUSSIAN_C, 
                                        cv2.THRESH_BINARY, 7, 8)

# Plot the images
images = [img, thresh_binary, adap_mean_2, adap_mean_2_inv,
          adap_mean_8, adap_gaussian_8]
fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(15, 15))

for ind, p in enumerate(images):
    ax = axs[ind % 2, ind // 2]
    ax.imshow(p, cmap='gray')
    ax.axis('off')

plt.show()

