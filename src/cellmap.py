"""
Functions for get the label image given the coords of the polygons' outer rings
"""
import pandas as pd
import json
import src.config as config
from scipy.sparse import coo_matrix, csr_matrix
import scipy.sparse as sp
import matplotlib.pyplot as plt
from skimage.measure import regionprops
from skimage.draw import line, polygon, circle, ellipse
import numpy as np
from PIL import Image, ImageDraw
import logging

logger = logging.getLogger()
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s:%(levelname)s:%(message)s"
)


def transformation(cgf):
    with open(cgf['manifest']) as f:
        settings = json.load(f)

    # bounding box in microns
    bbox = {'x0': settings['bbox_microns'][0],
            'x1': settings['bbox_microns'][2],
            'y0': settings['bbox_microns'][1],
            'y1': settings['bbox_microns'][3]}

    # image width and height in pixels
    img = {'width': settings['mosaic_width_pixels'],
           'height': settings['mosaic_height_pixels']}

    # Affine transformation: a set of coefficients a, b, c, d for transforming
    # a point of a form (x, y) into (a*x + b, c*y + d)
    a = img['width'] / (bbox['x1'] - bbox['x0'])
    b = -1 * img['width'] / (bbox['x1'] - bbox['x0']) * bbox['x0']
    c = img['height'] / (bbox['y1'] - bbox['y0'])
    d = -1 * img['height'] / (bbox['y1'] - bbox['y0']) * bbox['y0']

    tx = lambda x: a * x + b
    ty = lambda y: c * y + d
    return tx, ty, img, bbox


# def main():
#     cfg = config.DEFAULT
#     boundaries_file = "../dashboard/data/cellBoundaries.tsv"
#     boundaries = pd.read_csv(boundaries_file, sep='\t')
#     tx, ty, img, _ = transformation(cfg)
#
#     rr = []
#     cc = []
#     data = []
#     # boundaries = boundaries.head()
#     for index, row in boundaries.iterrows():
#         poly = row['coords']
#         poly = np.array(json.loads(poly))
#         x_px = tx(poly[:,0]).astype(np.uint32)
#         y_px = ty(poly[:,1]).astype(np.uint32)
#         _rr, _cc = polygon(x_px, y_px)
#         _data = np.ones(_cc.shape) * (index+1)
#         rr.extend(_rr.astype(np.uint32))
#         cc.extend(_cc.astype(np.uint32))
#         data.extend(_data)
#         logger.info('Doing cell %d out of %d' % (index, boundaries.shape[0]))
#
#
#     # rr = np.concatenate(rr)
#     # cc = np.concatenate(cc)
#     # data = np.concatenate(data)
#
#     coo = coo_matrix((data, (rr, cc)), dtype=np.uint32)
#     # coo = coo_matrix((data, (rr, cc)), shape=(img['height'], img['width']))
#     print('Label image has shape: ',  coo.shape)
#     print('ok')


def get_mask(boundaries, label):
    """
    returns the mask from a polygon's boundaries coords
    :return:
    """
    img = Image.new('L', (2000, 2000), 0)
    ImageDraw.Draw(img).polygon(boundaries, outline=label, fill=label)
    return np.array(img)


def get_label_image():
    cfg = config.DEFAULT
    boundaries_file = "../dashboard/data/cellBoundaries.tsv"
    outer_rings = pd.read_csv(boundaries_file, sep='\t')
    tx, ty, img, _ = transformation(cfg)

    coo_list = []
    coo_row_list = []
    coo_col_list = []
    coo_data_list = []
    df_list = []
    for index, row in outer_rings.iterrows():
        cell_id = row['cell_id']
        if index % 1000 == 0:
            logger.info('Doing cell id: %d from a total %d' % (cell_id, outer_rings.shape[0]))
        assert int(index + 1) == cell_id, 'Are these miss-alinged?'

        # 1. get the boundaries coords
        poly = row['coords']
        poly = np.array(json.loads(poly))

        # 2. Convert micron to pixel coords
        x_px = tx(poly[:,0]).astype(np.uint32)
        y_px = ty(poly[:,1]).astype(np.uint32)

        # 3. shift the coords
        offset_x = x_px.min()
        offset_y = y_px.min()
        x = x_px - offset_x
        y = y_px - offset_y

        # 4. fill now the polygon, label it. Make a mask
        boundaries = list(zip(x, y))
        mask = get_mask(boundaries, cell_id)

        # 5. Get area area and cell centroid
        props = regionprops(mask)[0]
        _df = pd.DataFrame({'cell_id': [cell_id],
                          'area': [props.area],
                          'centroid_x': [props.centroid[1] + offset_x],
                          'centroid_y': [props.centroid[0] + offset_y]})

        mask = coo_matrix(mask)

        # 5. make now a sparse matrix of the label image. Do not forget to shift back the coords
        _coo = coo_matrix((mask.data, (mask.row+offset_y, mask.col+offset_x)), shape=(img['height'], img['width']))
        coo_list.append(_coo)
        coo_row_list.append(_coo.row)
        coo_col_list.append(_coo.col)
        coo_data_list.append(_coo.data)
        df_list.append(_df)

    coo_row_stacked = np.hstack(coo_row_list)
    coo_col_stacked = np.hstack(coo_col_list)
    coo_data_stacked = np.hstack(coo_data_list)
    label_image = coo_matrix((coo_data_stacked, (coo_row_stacked, coo_col_stacked)), shape=(img['height'], img['width']))
    df = pd.concat(df_list, ignore_index=True)

    return label_image, df


if __name__ == "__main__":
    label_image, df = get_label_image()
    logger.info(label_image.shape)
