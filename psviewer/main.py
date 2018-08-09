import os
import argparse
from skimage import io
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons
import psviewer.terastitcher as ts
import psviewer.raw as raw


# Hold the state of the current FOV
c_idx = 0
z_idx = 0
y_idx = 0
x_idx = 0

cs = None
zs = None
ys = None
xs = None
fmt = None

zstep = None
pixel_width_um = None
clim = None
data_path = None
img = None

ax_img = None
fig, ax = plt.subplots(figsize=(8, 8))
plt.subplots_adjust(left=0.25, bottom=0.25)
ax_cmin = plt.axes([0.25, 0.1, 0.65, 0.03])
ax_cmax = plt.axes([0.25, 0.15, 0.65, 0.03])
ax_auto = plt.axes([0.8, 0.025, 0.1, 0.04])
ax_radio = plt.axes([0.025, 0.5, 0.15, 0.15])
slider_cmin = Slider(ax_cmin, 'Min Intensity', 0, 2**12-1, valinit=0, valstep=1, valfmt='%d')
slider_cmax = Slider(ax_cmax, 'Max Intensity', 0, 2**12-1, valinit=1000, valstep=1, valfmt='%d')
auto_btn = Button(ax_auto, 'Auto')
radio_btn = None


def adjust_cmin(event):
    global clim
    cmin = slider_cmin.val
    cmax = clim[1]
    clim = [cmin, cmax]
    if cmin < cmax:
        ax_img.set_clim(vmin=cmin, vmax=cmax)
    else:
        slider_cmin.set_val(cmax-1)
        ax_img.set_clim(vmin=cmax-1, vmax=cmax)


def adjust_cmax(event):
    global clim
    cmin = clim[0]
    cmax = slider_cmax.val
    clim = [cmin, cmax]
    if cmin < cmax:
        ax_img.set_clim(vmin=cmin, vmax=cmax)
    else:
        slider_cmax.set_val(cmin+1)
        ax_img.set_clim(vmin=cmin, vmax=cmin+1)


def auto_adjust_clim(event):
    global clim
    cmin = img.min()
    cmax = img.max()/2
    if cmin < cmax:
        clim = [cmin, cmax]
        slider_cmin.set_val(cmin)
        slider_cmax.set_val(cmax)


def change_channel(event):
    global c_idx
    new_c_idx = int(event[-1])
    if new_c_idx != c_idx:
        c_idx = new_c_idx
        plot_img()


slider_cmin.on_changed(adjust_cmin)
slider_cmax.on_changed(adjust_cmax)
auto_btn.on_clicked(auto_adjust_clim)


def plot_img():
    global ax_img, img
    c, z, y, x = cs[c_idx], zs[z_idx], ys[y_idx], xs[x_idx]
    img_path = ts.build_path(data_path, c, z, y, x, fmt)
    if fmt in ['tif', 'tiff']:
        img = io.imread(img_path)
    elif fmt == 'raw':
        img = raw.raw_imread(img_path)
    width = pixel_width_um * img.shape[0]
    x_um, y_um, z_um = x/10, y/10, z/10
    ax_img = ax.imshow(img, extent=[x_um, x_um+width, y_um+width, y_um], clim=clim, cmap='gray', interpolation='none')
    ax.set_title(f'{c}: z-slice ({z_idx+1}/{len(zs)}) {z_um} um')
    ax.set_xlabel(f'x-tile {x_idx+1} / {len(xs)}')
    ax.set_ylabel(f'y-tile {y_idx+1} / {len(ys)}')
    fig.canvas.draw()


def keypress(event):
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
        radio_btn.set_active(c_idx)
        plot_img()
    elif event.guiEvent.char in [str(i) for i in range(len(cs))]:  # catch numpad too
        c_idx = int(event.guiEvent.char)
        radio_btn.set_active(c_idx)
        plot_img()


def main():
    global cs, zs, ys, xs, fmt, pixel_width_um, clim, data_path, zstep, radio_btn
    parser = argparse.ArgumentParser(description='Pre-stitcher viewer')
    parser.add_argument('--data-path', type=str, help='input path to unstitcher terastitcher volume', default='.')
    parser.add_argument('--zstep', type=int, help='number of frames to move in z', default=100)
    parser.add_argument('--pixel-width', '-p', type=float, help='image pixel with in micron', default=1)
    parser.add_argument('--cmin', type=float, help='minimum of grayscale intensity display range', default=0)
    parser.add_argument('--cmax', type=float, help='maximum of grayscale intensity display range', default=1000)
    args = parser.parse_args()

    data_path = os.path.abspath(args.data_path)
    clim = [args.cmin, args.cmax]
    pixel_width_um = args.pixel_width
    zstep = args.zstep

    fmt = ts.get_format(data_path)
    cs, zs, ys, xs = ts.get_czyx(data_path)

    radio_btn = RadioButtons(ax_radio, cs, active=0)
    radio_btn.on_clicked(change_channel)
    fig.canvas.mpl_connect('key_press_event', keypress)
    plot_img()
    plt.show()

    # TODO: Allow changing the zstep size
    # TODO: cache images for faster rendering


if __name__ == '__main__':
    main()
