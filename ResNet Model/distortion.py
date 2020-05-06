#! /usr/bin/env python3

import cv2
import os
import pathlib
import random
import sys
import time
try:
    import cupy as np
    from cupy.fft import ifft2, fft2
    using_cupy = True
except:
    print('cupy not found, defaulting to numpy')
    import numpy as np
    from numpy.fft import ifft2, fft2
    using_cupy = False
import matplotlib.pyplot as plt
from aotools.turbulence.infinitephasescreen import PhaseScreenKolmogorov


PSK = None

def create_circular_mask(h, w, center=None, radius=None):
    # https://stackoverflow.com/questions/44865023/how-can-i-create-a-circular-mask-for-a-numpy-array

    if center is None: # use the middle of the image
        center = (w//2, h//2)
    if radius is None: # use the smallest distance between the center and image walls
        radius = min(center[0], center[1], w-center[0], h-center[1])

    Y, X = np.ogrid[:h, :w]
    dist_from_center = np.sqrt((X - center[0])**2 + (Y - center[1])**2)

    mask = dist_from_center <= radius
    return mask.astype(int)

def circular_aperture(img):
    orig_dim = img.shape

    pri_h = orig_dim[0] // 2 # 224 -> 112
    pri_w = orig_dim[1] // 2 # 224 -> 112

    sec_h = int(orig_dim[0] * 0.09375) # 224 -> 21
    sec_w = int(orig_dim[1] * 0.09375) # 224 -> 21
    sup_start = (pri_h-sec_h) // 2 # start boundary to super impose pupil2 mask into mask of size pupil1
    sup_end = sup_start + sec_h # end boundary to super impose pupil2 on pupil1

    pupil1 = create_circular_mask(pri_h, pri_w, radius=pri_h//2)
    pupil2_mask = create_circular_mask(sec_h, sec_w, radius=sec_h//2)
    pupil2 = np.zeros_like(pupil1)
    pupil2[sup_start:sup_end, sup_start:sup_end] = pupil2_mask

    pupil_mask = pupil1 - pupil2
    pupil_mask = np.pad(pupil_mask, ((orig_dim[0] - pri_h)// 2, ), 'constant')
    pupil_mask = pupil_mask * orig_dim[0] / np.sqrt(np.sum(pupil_mask ** 2))

    return pupil_mask

def atmospheric_distort(img, aperture_size=8, fried_param=0.164, outer_scale=100, random_seed=None, stencil_length_factor=32):
     '''Apply atmospheric distortion to input image
     Parameters:
         img (numpy.ndarray): Image to be distorted
         intensity (int): level of distortion 5 - minimal, 11 - medium, 20 - high
         aperture_size (int): size of capturing aperture in meters
         fried_param (float): size of atmospheric coherence length (Fried param) in meters
         outer_scale (int): It's in the AOTools Kolmogorov Phase Screen docs
         random_seed (int): seed to use to generate Kolmogorov Phase Screen
         stencil_length_factor (int): How much longer is stencil than desired phase

     Returns:
         numpy.ndarray: numpy array representing the distorted image 
     '''
     global PSK

     height, width = img.shape[0], img.shape[1]
     pxl_scale = aperture_size/height

     # error handling
     if height != width:
         return "Error! Height does not equal weight!"

     pupil_mask = circular_aperture(img)
     start = time.time()
     if not PSK:
         PSK = PhaseScreenKolmogorov(height, pxl_scale, fried_param, outer_scale, random_seed, stencil_length_factor)
     for i in range(random.randint(1, 75)):
         phase_screen = PSK.add_row()
     print('phase screen generated in {}s'.format(time.time() - start))
     if using_cupy:
         phase_screen = np.asarray(phase_screen)

     a = ifft2(pupil_mask * np.exp(1j * phase_screen))
     h = abs(abs(a) ** 2)
     img_slice = ifft2(fft2(h) * fft2(img[:,:,1])).real

     img_slice /= np.max(img_slice)
     img_slice *= 255
     return np.repeat(img_slice[:,:,np.newaxis], 3, axis=2)

def atmospheric_distort_image_file(filepath, output_directory, aperture_size=8, fried_param=0.164, outer_scale=100, random_seed=None, stencil_length_factor=32):
    def make_distorted_image_filename(filename):
        return '{}_{}_{}_{}_{}'.format(filename, aperture_size, fried_param, outer_scale, stencil_length_factor)

    img = cv2.imread(str(filepath))
    if using_cupy:
        img = np.asarray(img)
    img = atmospheric_distort(img, aperture_size, fried_param, outer_scale, random_seed, stencil_length_factor)
    if using_cupy:
        img = np.asnumpy(img)

    filename, ext = os.path.splitext(filepath.name)
    new_filename = make_distorted_image_filename(filename)
    new_filepath = output_directory/str(filepath).replace(filename, new_filename)

    new_filepath.parents[0].mkdir(parents=True, exist_ok=True)
    cv2.imwrite(str(new_filepath), img)
    return new_filepath

def atmospheric_distort_directory(directory_string, file_glob_matcher, output_directory, aperture_size=8, fried_param=0.164, outer_scale=100, random_seed=None, stencil_length_factor=32):
    #try:
    directory = pathlib.Path(directory_string)
    output_directory = pathlib.Path(output_directory)

    for path in directory.glob(file_glob_matcher):
            #try:
        print('attempting to distory image at: {}'.format(str(path)))
        start = time.time()
        new_file = atmospheric_distort_image_file(path, output_directory)
        end = time.time()
        print('distorted image created({}) at: {}'.format(end-start, str(new_file)))
            #except Exception as e:
            #    print('Unable to apply atmospheric distortion to file at: {} -- {}'.format(str(path), str(e)))
    #except Exception as e:
    #    print('Unable to process given directory: {}'.format(str(e)))
    #    raise e

if __name__ == '__main__':
    if len(sys.argv) != 4:
        print('Invalid command usage:')
        print('Proper usage (replace <> fields with values):')
        print('\tdistortion <directory relative path> <file matcher string> <output_directory>')

    try:
       directory_path = sys.argv[1]
       file_glob_matcher = sys.argv[2]
       output_directory = sys.argv[3]
       atmospheric_distort_directory(directory_path, file_glob_matcher, output_directory)
       print('Completed applying distortion to images')
    except Exception as e:
       print('Failed to apply distortion to images -- {}'.format(str(e)))
       raise e

