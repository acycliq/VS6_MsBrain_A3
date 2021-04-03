import datashader as ds
import numpy as np
from functools import partial
from datashader import transfer_functions as tf
from datashader.utils import export_image
import pandas as pd
import pyvips
import shutil
import os
import logging

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s:%(levelname)s:%(message)s"
)

logger = logging.getLogger()


def map_image_size(z):
    '''
    return the image size for each zoom level. Assumes that each map tile is 256x255
    :param z:
    :return:
    '''
    return 256 * 2 ** z


def tile_maker(z_depth, out_dir, img_path):
    # img_path = os.path.join(dir_path, 'demo_data', 'background_boundaries.tif')

    dim = map_image_size(z_depth)
    # remove the dir if it exists
    if os.path.exists(out_dir):
        shutil.rmtree(out_dir)

    # now make a fresh one
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    im = pyvips.Image.new_from_file(img_path, access='sequential')

    # The following two lines add an alpha component to rgb which allows for transparency.
    # Is this worth it? It adds quite a bit on the execution time, about x2 increase
    # im = im.colourspace('srgb')
    # im = im.addalpha()

    logger.info('Resizing image: %s' % img_path)
    factor = dim / max(im.width, im.height)
    im = im.resize(factor)
    logger.info('Done! Image is now %d by %d' % (im.width, im.height))
    pixel_dims = [im.width, im.height]

    # sanity check
    assert max(im.width, im.height) == dim, 'Something went wrong. Image isnt scaled up properly. ' \
                                            'It should be %d pixels in its longest side' % dim

    # im = im.gravity('south-west', dim, dim) # <---- Uncomment this if the origin is the bottomleft corner

    # now you can create a fresh one and populate it with tiles
    logger.info('Started doing the image tiles ')
    im.dzsave(out_dir, layout='google', suffix='.jpg', background=0, skip_blanks=0)
    logger.info('Done. Pyramid of tiles saved at: %s' % out_dir)

    return pixel_dims


def supertile(zoom_level, x_range, y_range):
    scene = ds.Canvas(plot_width=600, plot_height=600, x_range=x_range, y_range=y_range)
    aggregation = scene.points(data, 'x', 'y')
    image = tf.shade(aggregation, cmap=["#FF0000"], alpha=30)
    out = tf.spread(image, px=3, shape='circle', name="spread square")
    return out

np.random.seed(10)
mu = [0, 0]
scale = 100
cov = scale**2 * np.eye(2)
x, y = np.random.multivariate_normal(mu, cov, 50000).T

background = "red"
export = partial(export_image, background=None, export_path="exports")
x_range = [x.min(), x.max()]
y_range = [y.min(), y.max()]

data = pd.DataFrame({
    'x': x,
    'y': y
})
scene = ds.Canvas(plot_width=600, plot_height=600, x_range=x_range, y_range=y_range)
aggregation = scene.points(data, 'x', 'y')
image = tf.shade(aggregation, cmap=["#FF0000"], alpha=30)
out = tf.spread(image, px=3, shape='circle', name="spread square")
export(out, 'firery2b_600')
print('ok')