import os
import argparse
import multiprocessing
from skimage import io
from skimage.external import tifffile
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons
import psviewer.terastitcher as ts
import psviewer.raw as raw
from tqdm import tqdm


fig, ax = plt.subplots(figsize=(10, 10))


def stitch_img(input_path, c_idx, z_idx, shift_x, shift_y):
    fmt = ts.get_format(input_path)
    cs, zs, ys, xs = ts.get_czyx(input_path)

    c, z, y, x = cs[c_idx], zs[z_idx], ys[0], xs[0]
    fov_path = ts.build_path(input_path, c, z, y, x, fmt)
    if fmt in ['tif', 'tiff']:
        fov0 = io.imread(fov_path)
    elif fmt == 'raw':
        fov0 = raw.raw_imread(fov_path)
    fov_shape = fov0.shape

    dx = fov_shape[-1] - shift_x
    dy = fov_shape[1] - shift_y
    img_shape = (len(ys) * fov_shape[1], len(xs) * fov_shape[-1])

    img = np.zeros(img_shape, dtype=fov0.dtype)
    for y_idx in range(len(ys)):
        y_start = y_idx * dy
        y_stop = y_start + fov_shape[1]
        for x_idx in range(len(xs)):
            c, z, y, x = cs[c_idx], zs[z_idx], ys[y_idx], xs[x_idx]
            fov_path = ts.build_path(input_path, c, z, y, x, fmt)
            if fmt in ['tif', 'tiff']:
                fov = io.imread(fov_path)
            elif fmt == 'raw':
                fov = raw.raw_imread(fov_path)
            x_start = x_idx * dx
            x_stop = x_start + fov_shape[-1]
            prev_img = img[y_start:y_stop, x_start:x_stop]
            new_img = np.maximum(prev_img, fov)
            img[y_start:y_stop, x_start:x_stop] = new_img
    return img


def stitching_preview(input_path):
    c_idx, z_idx = 1, 1600
    shift_x, shift_y = 307, 307
    clim = [0, 600]
    img = stitch_img(input_path, c_idx, z_idx, shift_x, shift_y)
    plt.imshow(img, interpolation='none', animated=True, clim=clim)
    plt.show()


def _stitch_and_save(args):
    input_path, output_path, c_idx, z_idx, shift_x, shift_y, z = args
    filename = f'{z:06d}.tif'
    img = stitch_img(input_path, c_idx, z_idx, shift_x, shift_y)
    tifffile.imsave(os.path.join(output_path, filename), img, compress=1)


def stitch_color(input_path, output_path, c_idx, shift_x, shift_y, nb_workers=None):
    if nb_workers is None:
        nb_workers = multiprocessing.cpu_count()

    fmt = ts.get_format(input_path)
    cs, zs, ys, xs = ts.get_czyx(input_path)

    args_list = []
    for z_idx, z in enumerate(zs):
        args = (input_path, output_path, c_idx, z_idx, shift_y, shift_x, z)
        args_list.append(args)

    with multiprocessing.Pool(nb_workers) as pool:
            list(tqdm(pool.imap(_stitch_and_save, args_list), total=len(zs)))


def main():
    parser = argparse.ArgumentParser(description='Simple stitcher')
    parser.add_argument('--input-path', type=str, help='input path to unstitched images in terastitcher format', default='.')
    parser.add_argument('--output-path', type=str, help='output path for stitched images', default='stitched')
    args = parser.parse_args()

    input_path = os.path.abspath(args.input_path)
    output_path = os.path.abspath(args.output_path)
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    stitching_preview(input_path)

    c_idx = 1
    shift_x, shift_y = 307, 307
    stitch_color(input_path, output_path, c_idx, shift_x, shift_y)


if __name__ == '__main__':
    main()