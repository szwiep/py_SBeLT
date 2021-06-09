#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  4 21:49:03 2021

@author: smchartrand
"""


########################
from PIL import Image

# filepaths
fp_in = '../plots/iter{i}.png'
fp_out = '../plots/simulation_snapshot.gif'

length = 500
step = 1
images = []

for i in range(0, length, step):
    im = Image.open(fp_in.format(i=i))
    images.append(im)

images[0].save(fp_out, save_all=True, append_images=images[1:], optimize=True, duration=150, loop=0)

