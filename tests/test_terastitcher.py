import unittest
import os
import psviewer.terastitcher as ts


test_path = '/media/share2/SmartSPIM_Archive/2018_06_07/20180607_19_33_11_YG_PV-GPe_sample2_488LP75_561LP120_642LP75'
colors_true = ['Color_0', 'Color_1', 'Color_2']
xs_true = [77960, 109690, 141420, 173150]
ys_true = [2590, 34320, 66050, 97780, 129510, 161240]
zs_true = list(range(37210, 117191, 20))  # 2um step


class testFindFolders(unittest.TestCase):

    def test(self):
        folders = ts.find_folders(test_path)
        for color in colors_true:
            self.assertIn(color, folders)


class testFindImages(unittest.TestCase):

    def test(self):
        path = os.path.join(test_path, 'Color_0', f'{xs_true[0]:06d}', f'{xs_true[0]:06d}_{ys_true[0]:06d}')
        files = ts.find_images(path)
        self.assertEqual(len(files), len(zs_true))


class testFindColors(unittest.TestCase):

    def test_no_colors(self):
        folders = ts.find_folders(test_path)
        path = os.path.join(test_path, folders[0])
        colors = ts.find_colors(path)
        self.assertIsNone(colors)

    def test_3_colors(self):
        colors = ts.find_colors(test_path)
        self.assertEqual(colors, colors_true)


class testFindX(unittest.TestCase):

    def test_no_x(self):
        xs = ts.find_x(test_path)
        self.assertIsNone(xs)

    def test_x(self):
        folders = ts.find_folders(test_path)
        path = os.path.join(test_path, folders[0])
        xs = ts.find_x(path)
        self.assertEqual(xs, xs_true)


class testFindY(unittest.TestCase):

    def test_no_y(self):
        ys = ts.find_y(test_path)
        self.assertIsNone(ys)

    def test_y(self):
        folders = ts.find_folders(test_path)
        path_colors = os.path.join(test_path, folders[0])
        folders_x = ts.find_folders(path_colors)
        path = os.path.join(test_path, folders[0], folders_x[0])
        ys = ts.find_y(path)
        self.assertEqual(ys, ys_true)


class testFindZ(unittest.TestCase):

    def test_no_z(self):
        zs = ts.find_z(test_path)
        self.assertIsNone(zs)

    def test_z(self):
        path = os.path.join(test_path, 'Color_0', f'{xs_true[0]:06d}', f'{xs_true[0]:06d}_{ys_true[0]:06d}')
        zs = ts.find_z(path)
        self.assertEqual(zs, zs_true)


class testFindCZYX(unittest.TestCase):

    def test(self):
        cs, zs, ys, xs = ts.get_czyx(test_path)
        self.assertEqual(cs, colors_true)
        self.assertEqual(zs, zs_true)
        self.assertEqual(ys, ys_true)
        self.assertEqual(xs, xs_true)


class testBuildRelativePath(unittest.TestCase):

    def test(self):
        path = ts.build_relative_path(1, 2, 3, 4, 'raw')
        self.assertEqual(path, 'Color_1/000004/000004_000003/000002.raw')