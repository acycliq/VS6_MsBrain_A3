import pandas as pd
import numpy as np
import src.config as config
from src.get_data import is_inside_sm_parallel
import os
import logging

logger = logging.getLogger()
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s:%(levelname)s:%(message)s"
)


def get_spots(cfg):
    df_list = []
    for i in range(6):
        fName = os.path.join(cfg['target_dir'], 'geneData', 'geneData_%d.tsv' % i)
        df = pd.read_csv(fName, sep='\t')
        df_list.append(df)
    out = pd.concat(df_list, ignore_index=True)
    return out



if __name__ == "__main__":
    slice_id = 'MsBrain_Eg1_VS6_JH_V6_05-02-2021'
    region_id = 'region_0'
    cfg = config.get_config(slice_id=slice_id, region_id=region_id)
    # 1. read the rotated spots
    spots = get_spots(cfg)
    spots['inside_cell_key'] = np.zeros(spots.shape[0])
    logger.info(spots)

    # get the cell metadata
    cell_meta = pd.read_csv(cfg['cell_metadata'], index_col=0)

    # get the boundaries
    cellBoundaries = pd.read_csv(os.path.join(cfg['target_dir'], 'cellBoundaries', 'cellBoundaries.tsv'), sep='\t')
    cellProps = pd.read_csv(os.path.join(cfg['target_dir'], 'cell_props', 'cell_props.tsv'), sep='\t')

    # the boundaries appear to be strings. Convert them to a list of tuples
    _cell_boundaries = [eval(d) for d in cellBoundaries.cell_boundaries]
    xxx = []
    for i, x in enumerate(_cell_boundaries):
        _cell_boundaries[i] = [tuple(d) for d in x]
    cellBoundaries.cell_boundaries = _cell_boundaries

    # join them
    assert np.all(cellProps.cell_label == cellBoundaries.cell_label), 'cellBoundaries and cell_props are not aligned'
    cellBoundaries['cell_key'] = cellProps.cell_key
    cellBoundaries = cellBoundaries.set_index('cell_key')
    all_keys = set(cellBoundaries.index.values)
    # cellBoundaries = cellBoundaries.assign(cell_key=cellProps.cell_key.values)

    # 2. loop over the fovs
    fovs = np.unique(spots.fov.values)
    logger.info('start spot labelling')
    for fov in fovs:
        logger.info('Started fov: %d' % fov)
        # get the cell keys of the cells in the fov
        fov_cells_keys = cell_meta[cell_meta.fov == fov]

        # some of the cell keys may not be in the cell boundaries dataframe (probably because I looked on in z3 to capture the boundaries)
        missing = set(fov_cells_keys.index.values) - all_keys

        # remove the missing
        idx = set(fov_cells_keys.index.values) - missing

        # get the cells that are captured by this fov
        cells_in_fov = cellBoundaries.loc[list(idx)]

        spots_in_fov = spots[spots.fov == fov]
        points = spots_in_fov[['x', 'y']].values
        for cell_key, row in cells_in_fov.iterrows():
            outer_ring = row['cell_boundaries']
            if not outer_ring[-1] == outer_ring[0]:
                outer_ring.append(outer_ring[0])
            # find if a spots lies within the poly
            poly = np.array(outer_ring)
            mask = is_inside_sm_parallel(points, poly)

            # get the index of those spots that fall inside the poly
            spots_idx = spots_in_fov.index[mask]

            # backfill the inside_cell_key column with the corresponding key
            spots.loc[spots_idx, ['inside_cell_key']] = cell_key

        logger.info('spot labelling finished')
        print('ok')

    assert 1==0
    # 2. read the cell rotated polygons
    cells = pd.read_csv(r"E:\Neuroscience\dev\Python\VS6_MsBrain_A3\MsBrain_Eg1_VS6_JH_V6_05-02-2021\region_0\cellBoundaries\cellBoundaries_0.tsv", sep='\t')
    logger.info(cells)
    logger.info('Found %d cells' % cells.shape[0])

    # 3. get the roi
    if cfg['clip_poly']:
        logger.info('Found clipping poly. Keeping data inside %s' % cfg['clip_poly'])
        coords = cfg['clip_poly']
        if not coords[0] == coords[-1]:
            logger.info('Closing the polygon')
            coords.append(coords[0])

    # 3. loop over the polygon and label the spots
    for d, coords in enumerate(cells.coords):
        if not coords[0] == coords[-1]:
            coords = eval(coords)
            coords.append(coords[0])

        points = spots[['x', 'y']].values
        poly = np.array(coords)
        mask = is_inside_sm_parallel(points, poly)
        logger.info("Cell %d had %d spots" % (d, mask.sum()))
