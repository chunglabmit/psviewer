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
import pywt


# TODO: Allow shifts to be input and add a preview button
# TODO: Allow user to preview with different blending modes using radio buttons
# TODO: Implement wavelet-based blending mode and background subtraction


fig, ax = plt.subplots(figsize=(8, 8))


def _get_shape_and_dtype(input_path, c, z, y, x, fmt):
    fov_path = ts.build_path(input_path, c, z, y, x, fmt)
    if fmt in ['tif', 'tiff']:
        fov0 = io.imread(fov_path)
    elif fmt == 'raw':
        fov0 = raw.raw_imread(fov_path)
    return fov0.shape, fov0.dtype


def stitch_img_mip(input_path, c_idx, z_idx, shift_x, shift_y):
    fmt = ts.get_format(input_path)
    cs, zs, ys, xs = ts.get_czyx(input_path)
    c, z, y, x = cs[c_idx], zs[z_idx], ys[0], xs[0]
    fov_shape, fov_dtype = _get_shape_and_dtype(input_path, c, z, y, x, fmt)

    dx = fov_shape[-1] - shift_x
    dy = fov_shape[1] - shift_y
    img_shape = (len(ys) * fov_shape[1], len(xs) * fov_shape[-1])

    img = np.zeros(img_shape, dtype=fov_dtype)
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


def stitch_img_cosine(input_path, c_idx, z_idx, shift_x, shift_y):
    fmt = ts.get_format(input_path)
    cs, zs, ys, xs = ts.get_czyx(input_path)

    if z_idx == -1:
        z_idx = len(zs) // 2

    c, z, y, x = cs[c_idx], zs[z_idx], ys[0], xs[0]
    fov_shape, fov_dtype = _get_shape_and_dtype(input_path, c, z, y, x, fmt)

    # Open all images in this z-slice
    fovs = []
    for y_idx in range(len(ys)):
        fov_row = []
        for x_idx in range(len(xs)):
            c, z, y, x = cs[c_idx], zs[z_idx], ys[y_idx], xs[x_idx]
            fov_path = ts.build_path(input_path, c, z, y, x, fmt)
            if fmt in ['tif', 'tiff']:
                fov = io.imread(fov_path)
            elif fmt == 'raw':
                fov = raw.raw_imread(fov_path)
            fov_row.append(fov)
        fovs.append(fov_row)

    # Initialize the output image
    dx = fov_shape[-1] - shift_x
    dy = fov_shape[1] - shift_y
    img_shape = (len(ys) * fov_shape[0], len(xs) * fov_shape[-1])
    img = np.zeros(img_shape, dtype=fov_dtype)

    # Compute the blending functionpsviewer
    i_x = np.arange(shift_x)
    f_x = (1 + np.cos(i_x / (shift_x - 1) * np.pi)) / 2
    i_y = np.arange(shift_y)
    f_y = (1 + np.cos(i_y / (shift_y - 1) * np.pi)) / 2

    for y_idx, fov_row in enumerate(fovs):

        # copy the top row before blending the next row
        if y_idx != 0:
            row_top = np.copy(img[y_idx*dy:y_idx*dy+shift_y])

        img[y_idx*dy:y_idx*dy+fov_shape[0], :dx] = fov_row[0][:, :dx]

        for x_idx in range(1, len(fov_row)):
            fov_left = fov_row[x_idx-1]
            fov_right = fov_row[x_idx]

            # blend x intersections
            inter_left = fov_left[:, dx:]
            inter_right = fov_right[:, :shift_x]
            inter = f_x * inter_left + (1 - f_x) * inter_right

            # write the results
            img[y_idx*dy:y_idx*dy+fov_shape[0], x_idx*dx:x_idx*dx+shift_x] = inter
            img[y_idx*dy:y_idx*dy+fov_shape[0], x_idx*dx+shift_x:(x_idx+1)*dx] = fov_right[:, shift_x:dx]

        if y_idx == 0:  # No layer to blend above
            continue

        # blend y intersections
        row_btm = img[y_idx*dy:y_idx*dy+shift_y]
        inter = f_y[:, np.newaxis] * row_top + (1 - f_y[:, np.newaxis]) * row_btm

        # write the results
        img[y_idx*dy:y_idx*dy+shift_y] = inter

    return img


def wavelet_blend(left, right, shift_x, dx, f_x, wavelet, level):
    coeffs_left = pywt.swt2(left, wavelet, level)
    coeffs_right = pywt.swt2(right, wavelet, level)

    coeffs = []
    for l, (cleft, cright) in enumerate(zip(coeffs_left, coeffs_right)):
        approx_left, details_left = cleft
        approx_right, details_right = cright
        # blend approximation x intersections
        inter_approx_left = approx_left[:, dx:]
        inter_approx_right = approx_right[:, :shift_x]
        inter_approx = f_x * inter_approx_left + (1 - f_x) * inter_approx_right
        approx = approx_right
        approx[:, :shift_x] = inter_approx
        # blend detail x intersections
        details = []
        for dl, dr in zip(details_left, details_right):
            inter_dl = dl[:, dx:]
            inter_dr = dr[:, :shift_x]
            inter_detail = f_x * inter_dl + (1 - f_x) * inter_dr
            detail = dr
            detail[:, :shift_x] = inter_detail
            details.append(detail)
        coeffs.append((approx, tuple(details)))

    return pywt.iswt2(coeffs, wavelet)


def stitch_img_wavelet(input_path, c_idx, z_idx, shift_x, shift_y, wavelet='db5', level=4):
    fmt = ts.get_format(input_path)
    cs, zs, ys, xs = ts.get_czyx(input_path)
    c, z, y, x = cs[c_idx], zs[z_idx], ys[0], xs[0]
    fov_shape, fov_dtype = _get_shape_and_dtype(input_path, c, z, y, x, fmt)

    # Open all images in this z-slice
    fovs = []
    for y_idx in range(len(ys)):
        fov_row = []
        for x_idx in range(len(xs)):
            c, z, y, x = cs[c_idx], zs[z_idx], ys[y_idx], xs[x_idx]
            fov_path = ts.build_path(input_path, c, z, y, x, fmt)
            if fmt in ['tif', 'tiff']:
                fov = io.imread(fov_path)
            elif fmt == 'raw':
                fov = raw.raw_imread(fov_path)
            fov_row.append(fov)
        fovs.append(fov_row)

    # Initialize the output image
    dx = fov_shape[-1] - shift_x
    dy = fov_shape[1] - shift_y
    img_shape = (len(ys) * fov_shape[0], len(xs) * fov_shape[-1])
    img = np.zeros(img_shape)

    # Compute the blending function
    i_x = np.arange(shift_x)
    f_x = (1 + np.cos(i_x / (shift_x - 1) * np.pi)) / 2
    i_y = np.arange(shift_y)
    f_y = (1 + np.cos(i_y / (shift_y - 1) * np.pi)) / 2

    for y_idx, fov_row in enumerate(fovs):
        # copy the top row before blending the next row
        if y_idx != 0:
            row_top = np.copy(img[(y_idx-1)*dy:(y_idx-1)*dy+fov_shape[0]])

        img[y_idx * dy:y_idx * dy + fov_shape[0], :dx] = fov_row[0][:, :dx]

        for x_idx in range(1, len(fov_row)):
            fov_left = fov_row[x_idx - 1]
            fov_right = fov_row[x_idx]
            blended = wavelet_blend(fov_left, fov_right, shift_x, dx, f_x, wavelet, level)
            img[y_idx*dy:y_idx*dy+fov_shape[0], x_idx*dx:x_idx*dx+fov_shape[-1]] = blended

        if y_idx == 0:  # No layer to blend above
            continue

        # blend y intersections
        row_btm = np.copy(img[y_idx*dy:y_idx*dy+fov_shape[0]])

        for x_idx in range(len(fov_row)):
            fov_top = row_top[:, x_idx*dx:x_idx*dx+fov_shape[-1]]
            fov_btm = row_btm[:, x_idx*dx:x_idx*dx+fov_shape[-1]]
            blended = wavelet_blend(fov_top.T, fov_btm.T, shift_y, dy, f_y, wavelet, level)
            img[y_idx*dy:y_idx*dy+fov_shape[0], x_idx*dx:x_idx*dx+fov_shape[-1]] = blended.T

    return img


def stitching_preview(input_path, c_idx, z_idx, shift_x, shift_y, clim):
    img = stitch_img_cosine(input_path, c_idx, z_idx, shift_x, shift_y)

    plt.subplots_adjust(bottom=0.1)
    ax_img = ax.imshow(img,
                       clim=clim,
                       cmap='gray',
                       interpolation='none',
                       animated=True)
    ax_btn = plt.axes([0.66, 0.025, 0.1, 0.04])
    stitch_btn = Button(ax_btn, 'Stitch')
    stitch_btn.on_clicked(_close_and_stitch)
    plt.show()


def _stitch_and_save(args):
    input_path, output_path, c_idx, z_idx, shift_x, shift_y, z = args
    filename = f'{z:06d}.tif'
    img = stitch_img_cosine(input_path, c_idx, z_idx, shift_x, shift_y)
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


stitch_flag = False


def _close_and_stitch(event):
    global stitch_flag
    plt.close()
    stitch_flag = True


def main():
    parser = argparse.ArgumentParser(description='Simple stitcher')
    parser.add_argument('--input-path', type=str, help='input path to unstitched images in terastitcher format', default='.')
    parser.add_argument('--output-path', type=str, help='output path for stitched images', default='stitched')
    parser.add_argument('--c-idx', '-c', type=int, help='color index to preview stitching', default=0)
    parser.add_argument('--z-idx', '-z', type=int, help='z index to preview stitching', default=-1)
    parser.add_argument('--shift-x', type=int, help='number of pixel overlap between tiles in x', default=307)
    parser.add_argument('--shift-y', type=int, help='number of pixel overlap between tiles in y', default=307)
    parser.add_argument('--cmin', type=float, help='minimum intensity for display range', default=0)
    parser.add_argument('--cmax', type=float, help='maximum intensity for display range', default=600)
    args = parser.parse_args()

    input_path = os.path.abspath(args.input_path)
    output_path = os.path.abspath(args.output_path)
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    c_idx = args.c_idx
    z_idx = args.z_idx
    shift_x, shift_y = args.shift_x, args.shift_y
    cmin = args.cmin
    cmax = args.cmax
    clim = [cmin, cmax]

    stitching_preview(input_path, c_idx, z_idx, shift_x, shift_y, clim)

    if stitch_flag is True:
        stitch_color(input_path, output_path, c_idx, shift_x, shift_y)
    else:
        print('Viewer closed, stitching aborted')


if __name__ == '__main__':
    main()