"""
Convenience functions only to collect the data needed for the draft viewer
It reads the vizgen raw data and creates geneData, cellData and cellBoundaries needed
to visualise the data. It DOES NOT DO ANY CELLTYPING
"""

import os
import numpy as np
import pandas as pd
import src.config as config
import logging
from geopandas import GeoSeries
from shapely.geometry import Point, Polygon
import numba
from src.utils import transformation
from src.utils import splitter_mb, save_df
from src.cellBorders import cell_boundaries_px, cell_boundaries_px_par
import src.cellBorders as cellBorders
from src.utils import rotate_data

from contextlib import closing # this will correctly close the request
import io
# import dropbox


logger = logging.getLogger()
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s:%(levelname)s:%(message)s"
)


# def clip_data(df, cfg):
#     """
#     Keeps spots inside the area defined by predetermined polygon only
#     :param df:
#     :param cfg:
#     :return:
#     """
#     if cfg['clip_poly']:
#         logger.info('Found clipping poly. Keeping data inside %s' % cfg['clip_poly'])
#         coords = cfg['clip_poly']
#         # logger.info('Rotating the clipping polygon by %d degrees' % cfg['rotation'][0])
#         # coords_rot = rotate_data(np.array(coords), cfg).astype(np.int32)
#         # logger.info('Coords of the rotated polygon are %s ' % coords_rot.tolist())
#         # poly = Polygon(coords_rot)
#         poly = Polygon(coords)
#
#         s = GeoSeries(map(Point, zip(df.x.values, df.y.values)))
#         mask = s.within(poly).values
#         return df[mask]
#     else:
#         return df


def clip_data(df, cfg):
    """
    Keeps spots inside the area defined by predetermined polygon only
    :param df:
    :param cfg:
    :return:
    """
    if cfg['clip_poly']:
        logger.info('Found clipping poly. Keeping data inside %s' % cfg['clip_poly'])
        coords = cfg['clip_poly']
        if not coords[0] == coords[-1]:
            logger.info('Closing the polygon')
            coords.append(coords[0])

        points = df[['x', 'y']].values
        poly = np.array(coords)
        mask = is_inside_sm_parallel(points, poly)
        logger.info('Collected %d points inside ROI' % mask.sum())
        return df[mask]
    else:
        return df


@numba.jit(nopython=True)
def is_inside_sm(polygon, point):
    # From https://github.com/sasamil/PointInPolygon_Py/blob/master/pointInside.py
    # and
    # https://github.com/sasamil/PointInPolygon_Py/blob/master/pointInside.py
    length = len(polygon)-1
    dy2 = point[1] - polygon[0][1]
    intersections = 0
    ii = 0
    jj = 1

    while ii<length:
        dy  = dy2
        dy2 = point[1] - polygon[jj][1]

        # consider only lines which are not completely above/bellow/right from the point
        if dy*dy2 <= 0.0 and (point[0] >= polygon[ii][0] or point[0] >= polygon[jj][0]):

            # non-horizontal line
            if dy<0 or dy2<0:
                F = dy*(polygon[jj][0] - polygon[ii][0])/(dy-dy2) + polygon[ii][0]

                if point[0] > F: # if line is left from the point - the ray moving towards left, will intersect it
                    intersections += 1
                elif point[0] == F: # point on line
                    return 2

            # point on upper peak (dy2=dx2=0) or horizontal line (dy=dy2=0 and dx*dx2<=0)
            elif dy2==0 and (point[0]==polygon[jj][0] or (dy==0 and (point[0]-polygon[ii][0])*(point[0]-polygon[jj][0])<=0)):
                return 2

        ii = jj
        jj += 1

    #print 'intersections =', intersections
    return intersections & 1

@numba.njit(parallel=True)
def is_inside_sm_parallel(points, polygon):
    ln = len(points)
    D = np.empty(ln, dtype=numba.boolean)
    for i in numba.prange(ln):
        D[i] = is_inside_sm(polygon,points[i])
    return D


def get_gene_data(cfg):
    # Transformation from micron to pixel coords
    tx, ty, _, _ = transformation(cfg)

    # read the spots from the raw files
    spots_path = cfg['detected_transcripts']
    chunks = pd.read_csv(spots_path, chunksize=100000)
    data = pd.concat(chunks)

    # data_z3 = data[data.global_z == 3]
    data_z3 = data
    [unique_gene_names, gene_id, counts_per_gene] = np.unique(data_z3.gene.values, return_inverse=True,
                                                              return_counts=True)

    data_z3 = data_z3.assign(global_x_px=tx(data_z3.global_x.values).astype(np.int32))
    data_z3 = data_z3.assign(global_y_px=ty(data_z3.global_y.values).astype(np.int32))
    data_z3['Gene_id'] = gene_id
    data_z3 = data_z3.sort_values(by=['global_x_px', 'global_y_px'])

    data_z3 = data_z3[['gene', 'global_x_px', 'global_y_px', 'global_z', 'Gene_id']].rename(
        columns={'gene': "Gene", 'global_x_px': "x", 'global_y_px': "y", 'global_z': 'z_stack'})

    data_z3['neighbour'] = np.ones(len(gene_id)).astype(np.int32)
    data_z3['neighbour_array'] = [[1] for i in range(len(gene_id))]
    data_z3['neighbour_prob'] = [[1.0] for i in range(len(gene_id))]

    # data_z3 = data_z3.drop([['Unnamed: 0', 'Unnamed: 0.1']], axis=1)

    logger.info('Data size: %d, %d' % (data_z3.shape[0], data_z3.shape[1]))
    data_z3 = data_z3.sort_values(by=['x', 'y'])
    # splitter_mb(data_z3, 'geneData', 99)
    return data_z3


def run(slice_id, region_id):
    cfg = config.get_config(slice_id=slice_id, region_id=region_id)
    out_path = os.path.join('D:\\rotated_dapi_map_tiles', slice_id, region_id)

    # 1. First fetch the data
    geneData = get_gene_data(cfg)

    # 2. Keep only spots inside the ROI
    logger.info('Keeping only the ROI spots')
    geneData = clip_data(geneData, cfg)

    # 3. Rotate the spots
    logger.info('Rotating the transcript data by %d degrees' % cfg['rotation'][0])
    rot = rotate_data(geneData[['x', 'y']].values.copy(), cfg)
    geneData.x = rot[:, 0].astype(np.int32)
    geneData.y = rot[:, 1].astype(np.int32)

    # 4. Save the spots to the disk
    # splitter_mb(geneData, os.path.join(out_path, 'geneData'), 99)
    save_df(geneData, os.path.join(out_path, 'geneData'))
    logger.info('Gene data saved at: %s' % os.path.join(out_path, 'geneData'))

    # px_boundaries = cell_boundaries_px_par(cfg)
    # boundaries = px_boundaries[['cell_label', 'cell_boundaries']]
    # cellBorders.write_tsv(boundaries, os.path.join(out_path))
    # logger.info('cell data saved at: %s' % os.path.join(out_path))


if __name__ == "__main__":
    slice_ids = [
        "MsBrain_Eg1_VS6_JH_V6_05-02-2021",
        "MsBrain_Eg2_VS6_V11_JH_05-02-2021",
        "MsBrain_Eg3_VS6_JH_V6_05-01-2021",
        "MsBrain_EG4_VS6library_V6_LH_04-14-21",
        "MsBrain_Eg5_VS6_JH_V6_05-16-2021"
        ]
    region_ids = ['region_0', 'region_1']

    for slice_id in slice_ids:
        for region_id in region_ids:
            logger.info("\n Started slice %s, region %s" % (slice_id, region_id))
            try:
                run(slice_id, region_id)
            except KeyError as e:
                logger.info('KeyError %s' % str(e))
            except FileNotFoundError as e:
                logger.info('FileNotFoundError %s' % str(e))

    logger.info('Done!')
