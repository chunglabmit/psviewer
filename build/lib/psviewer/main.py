import os
import argparse
from skimage import io
import numpy as np
import matplotlib.pyplot as plt
import psviewer.terastitcher as ts
import psviewer.raw as raw


zstep = 200

c_idx = 0
z_idx = 0
y_idx = 0
x_idx = 0


def plot_img():
    c, z, y, x = cs[c_idx], zs[z_idx], ys[y_idx], xs[x_idx]
    img_path = ts.build_path(data_path, c, z, y, x, fmt)
    if fmt in ['tif', 'tiff']:
        img = io.imread(img_path)
    elif fmt == 'raw':
        img = raw.raw_imread(img_path)
    width = pixel_width_um * img.shape[0]
    x_um, y_um, z_um = x/10, y/10, z/10
    ax.imshow(img, extent=[x_um, x_um+width, y_um+width, y_um], clim=clim, cmap='gray', interpolation='none')
    ax.set_title(f'{c}: z-slice ({z_idx+1}/{len(zs)}) {z_um} um')
    ax.set_xlabel(f'x-tile {x_idx+1} / {len(xs)}')
    ax.set_ylabel(f'y-tile {y_idx+1} / {len(ys)}')
    fig.canvas.draw()


def press(event):
    global c_idx, z_idx, y_idx, x_idx
    if event.key == ',':
        if z_idx >= zstep:
            z_idx -= zstep
        else:
            z_idx = 0
        plot_img()
    elif event.key == '.':
        if z_idx < len(zs)-zstep:
            z_idx += zstep
        else:
            z_idx = len(zs)-1
        plot_img()
    elif event.key == 'up':
        if y_idx >= 1:
            y_idx -= 1
        plot_img()
    elif event.key == 'down':
        if y_idx < len(ys)-1:
            y_idx += 1
        plot_img()
    elif event.key == 'left':
        if x_idx >= 1:
            x_idx -= 1
        plot_img()
    elif event.key == 'right':
        if x_idx < len(xs)-1:
            x_idx += 1
        plot_img()
    elif event.key in [str(i) for i in range(len(cs))]:
        c_idx = int(event.key)
        plot_img()
    # print(f'Color_{c_idx}: z = {z_idx}, y = {y_idx}, x = {x_idx}')


def main():
    parser = argparse.ArgumentParser(description='Pre-stitcher viewer')
    parser.add_argument('data_path', type=str, help='input path to unstitcher terastitcher volume')
    parser.add_argument('--pixel-width', '-p', type=float, help='image pixel with in micron', default=1.0)
    parser.add_argument('--cmin', type=float, help='minimum of grayscale intensity display range')
    parser.add_argument('--cmax', type=float, help='maximum of grayscale intensity display range')
    args = parser.parse_args()

    data_path = os.path.abspath(args['data_path'])
    clim = [args['cmin'], args['cmax']]
    pixel_width_um = args['pixel-width']

    print(pixel_width_um)



# fmt = ts.get_format(data_path)
# cs, zs, ys, xs = ts.get_czyx(data_path)
#
# fig, ax = plt.subplots(figsize=(8, 8))
# fig.canvas.mpl_connect('key_press_event', press)
# plot_img()
# plt.show()

# TODO: Allow changing the zstep size
# TODO: cache images for faster rendering

if __name__ == '__main__':
    main()
    # data_path = '/media/share2/SmartSPIM_Archive/2018_06_07/20180607_19_33_11_YG_PV-GPe_sample2_488LP75_561LP120_642LP75'
    # clim = [0, 800]
    # pixel_width_um = 1.6

