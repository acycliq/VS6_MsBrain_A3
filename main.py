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
import json
from scipy.sparse import coo_matrix, save_npz, load_npz
from src.cellBorders import cell_boundaries_px
from src.cellmap import get_label_image, get_label_image_par
from src.utils import transformation

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
    data_z3 = data_z3.assign(global_x_px=tx(data_z3.global_x.values))
    data_z3 = data_z3.assign(global_y_px=ty(data_z3.global_y.values))
    data_z3 = data_z3.sort_values(by=['global_x_px', 'global_y_px'])
    data_z3 = data_z3.drop(['Unnamed: 0'], axis=1)
    return data_z3


def label_spots(label_image, spots: np.array) -> np.array:
    spots = spots.astype(np.int)
    m = label_image[spots[:, 1], spots[:, 0]]
    out = np.asarray(m)
    return out[0]


def run(cfg):
    CACHED_BOUNDARIES = False

    if CACHED_BOUNDARIES:
        boundaries = pd.read_csv('cell_boundaries_px.csv')
        boundaries['cell_boundaries'] = boundaries.cell_boundaries.apply(json.loads)
    else:
        px_boundaries = cell_boundaries_px(cfg)
        px_boundaries.to_csv('cell_boundaries_px.csv')

    # label_image, cell_props = get_label_image(px_boundaries, cfg)
    label_image, cell_props = get_label_image(px_boundaries, cfg)
    save_npz('coo_label_image.npz', coo_matrix(label_image))
    spots_df = read_spots(cfg)

    spots = spots_df[['global_x_px', 'global_y_px']].values
    labels = label_spots(label_image.tocsc(), spots)
    spots_df = spots_df.assign(cell_label=labels)

    # save now the results
    spots_df.to_csv('transcripts.csv', index=False)
    cell_props.to_csv('cell_props.csv', index=False)


if __name__=="__main__":
    cfg = config.DEFAULT
    run(cfg)
    logger.info('Done!')




#
# barcode_id = data.barcode_id.values
# nG = len(np.unique(barcode_id))
# cell_label = data.cell_label.values
# nC = len(np.unique(cell_label))
# group_idx = np.array([cell_label, barcode_id])
# df = pd.DataFrame(npg.aggregate(group_idx, barcode_id, size=[355710, 315], func='count'))
# df.head()
#
# # Get rid of the background
# df = df.iloc[1:, :]
# df