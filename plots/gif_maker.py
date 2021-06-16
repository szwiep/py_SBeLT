#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  4 21:49:03 2021

@author: smchartrand
"""


########################
from PIL import Image

# filepaths
in_dir = '../plots/shelf_test/'
in_filename = 'iter{i}.png'

out_dir = '../plots/shelf_test/'
out_filename = 'simulation_snapshot.gif'

fp_in = in_dir + in_filename
fp_out = out_dir + out_filename

start = 0
stop = 99
step = 1
images = []

for i in range(start, stop, step):
    im = Image.open(fp_in.format(i=i))
    images.append(im)

images[0].save(fp_out, save_all=True, append_images=images[1:], optimize=True, duration=150, loop=0)

