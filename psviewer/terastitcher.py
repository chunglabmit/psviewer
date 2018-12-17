"""This module searches for image according to the terastitcher naming conventions

Terastitcher input volumes follow a folder structure corresponding to:
    Color_C > XXXXXX > XXXXXX_YYYYYY > ZZZZZZ.tif
This module provides tools for parsing these volumes
"""
import os


MAX_CHANNELS = 10


def find_folders(path):
    return next(os.walk(path))[1]


def _check_empty_then_sort(coll):
    if not coll:
        return None
    coll.sort()
    return coll


def find_images(path):
    files = []
    for file in os.listdir(path):
        if file.endswith(('.tif', '.tiff', '.raw')):
            files.append(file)
    return _check_empty_then_sort(files)


def find_colors(path):
    folders = find_folders(path)
    colors = []
    for folder in folders:
        if folder in [f'Color_{i}' for i in range(MAX_CHANNELS)]:
            colors.append(folder)
    if len(colors) == 0:  # Add support for new naming convention
        for folder in folders:
            if folder in [f'Ex_{i}_Em_{i}' for i in range(MAX_CHANNELS)]:
                colors.append(folder)
    return _check_empty_then_sort(colors)


def find_x(path):
    folders = find_folders(path)
    xs = []
    for folder in folders:
        if len(folder) == 6:  # xxxxxx
            xs.append(int(folder))
    return _check_empty_then_sort(xs)


def find_y(path):
    folders = find_folders(path)
    ys = []
    for folder in folders:
        if len(folder) == 13:  # xxxxxx_yyyyyy
            ys.append(int(folder[-6:]))
    return _check_empty_then_sort(ys)


def find_z(path):
    filenames = find_images(path)
    if not filenames:  # skip if not images found
        return None
    zs = []
    for filename in filenames:
        basename = filename.split('.')[0]
        if len(basename) == 6:  # zzzzzz
            zs.append(int(basename))
    return _check_empty_then_sort(zs)


def get_czyx(path):
    cs = find_colors(path)
    tmp_path = os.path.join(path, cs[0])
    xs = find_x(tmp_path)
    tmp_path = os.path.join(tmp_path, f'{xs[0]:06d}')
    ys = find_y(tmp_path)
    tmp_path = os.path.join(tmp_path, f'{xs[0]:06d}_{ys[0]:06d}')
    zs = find_z(tmp_path)
    return cs, zs, ys, xs


def get_format(path):
    cs = find_colors(path)
    tmp_path = os.path.join(path, cs[0])
    xs = find_x(tmp_path)
    tmp_path = os.path.join(tmp_path, f'{xs[0]:06d}')
    ys = find_y(tmp_path)
    tmp_path = os.path.join(tmp_path, f'{xs[0]:06d}_{ys[0]:06d}')
    filenames = find_images(tmp_path)
    return filenames[0].split('.')[-1]


def build_relative_path(c, z, y, x, fmt):
    return os.path.join(c, f'{x:06d}', f'{x:06d}_{y:06d}', f'{z:06d}.{fmt}')


def build_path(basepath, c, z, y, x, fmt):
    return os.path.join(basepath, build_relative_path(c, z, y, x, fmt))
