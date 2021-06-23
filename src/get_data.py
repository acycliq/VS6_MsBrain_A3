"""
Convenience functions only to collect the data needed for the draft viewer
It reads the vizgen raw data and creates geneData, cellData and cellBoundaries needed
to visualise the data. It DOES NOT DO ANY CELLTYPING
"""
import sys
import os
import numpy as np
import pandas as pd
import src.config as config
import logging
from geopandas import GeoSeries
from shapely.geometry import Point, Polygon
import credentials
from src.utils import transformation
from src.utils import splitter_mb, dropbox_streamer
from src.cellBorders import cell_boundaries_px
import src.cellBorders as cellBorders

from contextlib import closing # this will correctly close the request
import io
# import dropbox


logger = logging.getLogger()
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s:%(levelname)s:%(message)s"
)


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
        poly = Polygon(coords)

        s = GeoSeries(map(Point, zip(df.x.values, df.y.values)))
        mask = s.within(poly).values
        return df[mask]
    else:
        return df


def get_gene_data(cfg):
    # Transormation from micron to pixel coords
    tx, ty, _, _ = transformation(cfg)

    # read the spots from the raw files
    spots_path = cfg['detected_transcripts']
    chunks = pd.read_csv(spots_path, chunksize=100000)
    data = pd.concat(chunks)
    # data = dropbox_streamer(spots_path)
    data_z3 = data[data.global_z == 3]
    [unique_gene_names, gene_id, counts_per_gene] = np.unique(data_z3.gene.values, return_inverse=True,
                                                              return_counts=True)

    data_z3 = data_z3.assign(global_x_px=tx(data_z3.global_x.values).astype(np.int32))
    data_z3 = data_z3.assign(global_y_px=ty(data_z3.global_y.values).astype(np.int32))
    data_z3 = data_z3.sort_values(by=['global_x_px', 'global_y_px'])

    data_z3 = data_z3[['gene', 'global_x_px', 'global_y_px']].rename(
        columns={'gene': "Gene", 'global_x_px': "x", 'global_y_px': "y"})

    data_z3['Gene_id'] = gene_id
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
    out_path = os.path.join(slice_id, region_id)

    geneData = get_gene_data(cfg)
    geneData = clip_data(geneData, cfg)
    splitter_mb(geneData, os.path.join(out_path, 'geneData'), 99)
    logger.info('Gene data saved at: %s' % os.path.join(out_path, 'geneData'))

    # px_boundaries = cell_boundaries_px(cfg)
    # boundaries = px_boundaries[['cell_label', 'cell_boundaries']]
    # cellBorders.write_tsv(boundaries, os.path.join(out_path))
    # logger.info('cell data saved at: %s' % os.path.join(out_path))


if __name__ == "__main__":
    slice_ids = [
        "MsBrain_Eg1_VS6_JH_V6_05-02-2021",
        "MsBrain_Eg2_VS6_V11_JH_05-02-2021",
        "MsBrain_Eg3_VS6_JH_V6_05-01-2021",  # missing file region_1\\images\\manifest.json
        "MsBrain_EG4_VS6library_V6_LH_04-14-21", # missing file region_1\\images\\manifest.json
        "MsBrain_Eg5_VS6_JH_V6_05-16-2021"
        ]
    region_ids = ['region_0', 'region_1']

    for slice_id in slice_ids:
        logger.info('Started slice %s' % slice_id)
        for region_id in region_ids:
            logger.info("Started region %s" % region_id)
            run(slice_id, region_id)

    logger.info('Done!')
