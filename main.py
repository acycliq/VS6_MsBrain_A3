"""
Main script to:
 1. get the cell boundaries
 2. get label_image
 3. label the spots
"""
import numpy as np
import pandas as pd
import src.config as config
import logging
from src.cellBorders import cell_boundaries
from src.cellmap import get_label_image, transformation

logger = logging.getLogger()
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s:%(levelname)s:%(message)s"
)


def read_spots(cfg):
    # Transormation from micron to pixel coords
    tx, ty, _, _ = transformation(cfg)
    # read the spots from the raw files
    spots_path = cfg['detected_transcripts']
    chunks = pd.read_csv(spots_path, chunksize=100000)
    data = pd.concat(chunks)
    data_z3 = data[data.global_z == 3]

    logger.info('Data size: %d, %d' % (data_z3.shape[0], data_z3.shape[1]))
    data_z3 = data_z3.assign(global_x_px = tx(data_z3.global_x.values))
    data_z3 = data_z3.assign(global_y_px= tx(data_z3.global_y.values))
    data_z3 = data_z3.sort_values(by=['global_x_px', 'global_y_px'])
    data_z3 = data_z3.drop(['Unnamed: 0'], axis=1)
    return data_z3


def label_spots(label_image, spots: np.array) -> np.array:
    spots = spots.astype(np.int)
    m = label_image[spots[:, 1], spots[:, 0]]
    out = np.asarray(m)
    return out[0]


def run(cfg):
    boundaries = cell_boundaries()
    label_image, cell_props = get_label_image(boundaries.head(10000), cfg)
    spots_df = read_spots(cfg)

    spots = spots_df[['global_x_px', 'global_y_px']].values
    spot_labels = label_spots(label_image.tocsc(), spots)
    spots_df = spots_df.assign(spots_label=spot_labels)

    print('OK')


if __name__=="__main__":
    cfg = config.DEFAULT
    run(cfg)