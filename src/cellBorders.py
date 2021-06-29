"""
Functions for get the coordinates of the outer rings describing the cell boundaries
"""
import os
import h5py
from os import listdir
from os.path import isfile, join
import numpy as np
import pandas as pd
from natsort import natsorted
from multiprocessing import Pool
from functools import partial
from itertools import chain
import src.config as config
from src.utils import splitter_mb, transformation
from src.utils import rotate_data
import logging

logger = logging.getLogger()
# logging.basicConfig(
#     level=logging.INFO,
#     format="%(asctime)s:%(levelname)s:%(message)s"
# )


def write_boundaries_tsv(boundaries, target_dir):
    N = len(boundaries)
    dummy_list = [[] for d in boundaries.cell_label]
    dummy_prob = [[1.0] for d in boundaries.cell_label]
    dummy_name = [['Generic'] for d in boundaries.cell_label]
    df = pd.DataFrame({'Cell_Num': boundaries.cell_label.values,
                       'X': np.zeros(N).astype(np.int),
                       'Y': np.zeros(N).astype(np.int),
                       'Genenames': dummy_list,
                       'CellGeneCount': dummy_list,
                       'ClassName': dummy_name,
                       'Prob': dummy_prob
                       })
    # df.to_csv(os.path.join(target_dir, 'cellData.tsv'), sep='\t', index=False)
    splitter_mb(df, os.path.join(target_dir, 'cellData'), 99)
    logger.info('cellData saved at: %s' % os.path.join(target_dir, 'cellData'))


def write_celldata_tsv(boundaries, target_dir):
    df = pd.DataFrame({'cell_id': boundaries.cell_label.values,
                       'coords': boundaries.cell_boundaries})
    df.to_csv(os.path.join(target_dir, 'cellBoundaries.tsv'), sep='\t', index=False)
    splitter_mb(df, os.path.join(target_dir, 'cellBoundaries'), 99)
    logger.info('cellBoundaries saved at: %s' % os.path.join(target_dir, 'cellBoundaries'))


def write_tsv(boundaries, out_dir):
    write_boundaries_tsv(boundaries, out_dir)
    write_celldata_tsv(boundaries, out_dir)


def micron_to_pixel(temp, cfg):
    tx, ty, _, _ = transformation(cfg)
    x = tx(temp[:, 0]).astype(np.int)
    y = ty(temp[:, 1]).astype(np.int)
    x_max = x.max()
    y_max = y.max()
    return np.array(list(zip(x, y))).tolist(), x_max, y_max


def get_boundaries(fileName: str, cfg, zIndex:int=3):
    """
    Helper function that reads a cached hdf5 file and gets the cell boundaries for that specific fov only
    :param fileName:
    :param zIndex:
    :return:
    """
    # target_hdf5 = os.path.join(config.ROOT_DIR, 'cell_boundaries', fileName)
    target_hdf5 = os.path.join(cfg['cell_boundaries_dir'], fileName)
    with h5py.File(target_hdf5, 'r') as data:
        cell_boundaries = []
        cell_keys = []

        for cell_key in data['featuredata'].keys():
            for p in data['featuredata'][cell_key]['zIndex_%i' % zIndex]:
                temp = np.array(data['featuredata'][cell_key]['zIndex_%i' % zIndex][p]['coordinates'][0])
                temp = outline_min(temp)
                temp, x_max, y_max = micron_to_pixel(temp, cfg)
                cell_boundaries.append(temp)
                cell_keys.append(cell_key)
    return cell_boundaries, cell_keys


def outline_min(arr: np.array) -> np.array:
    # minifies the boundaries array by casting the coords as int and then
    # removing the duplicates
    arr = arr.astype(np.int)
    indexes = np.unique(arr, axis=0, return_index=True)[1]
    arr_min = [arr[index].tolist() for index in sorted(indexes)]
    arr_min.append(arr_min[0])  # close the outline by appending the first point at the end
    return np.array(arr_min)


def cell_boundaries_px(cfg):
    """
    :return: returns a list of lists. Each sublist is the list of the x,y coords describing the cell boundaries
    """
    # hdf5_dir = os.path.join(config.ROOT_DIR, 'cell_boundaries')
    hdf5_dir = cfg['cell_boundaries_dir']
    # metadata_csv = config.DEFAULT['cell_metadata']
    # cellMetadata = pd.read_csv(metadata_csv).rename(columns={'Unnamed: 0': 'uid'})

    hfd5_files = [f for f in listdir(hdf5_dir) if isfile(join(hdf5_dir, f)) and os.path.splitext(f)[1] == '.hdf5']
    boundaries = []
    keys = []
    for hdf5_file in natsorted(hfd5_files):
        logger.info('rotating all cell boundaries in %s by %d degrees' % (hdf5_file, cfg['rotation'][0]))
        fov_cells, fov_cell_keys = get_boundaries(hdf5_file, cfg)
        for fov_cell, fov_cell_key in zip(fov_cells, fov_cell_keys):
            fov_cell = rotate_data(fov_cell, cfg)
            boundaries.append(fov_cell)
            keys.append(fov_cell_key)
    out = pd.DataFrame({'cell_key': keys,
                        'cell_label': np.arange(1, len(boundaries) + 1).astype(np.int),
                        'cell_boundaries': boundaries})
    # out.to_csv('cell_boundaries.csv', index=False)
    return out


def cell_boundaries_px_par(cfg):
    """
    :return: returns a list of lists. Each sublist is the list of the x,y coords describing the cell boundaries
    """
    hdf5_dir = cfg['cell_boundaries_dir']
    hfd5_files = [f for f in listdir(hdf5_dir) if isfile(join(hdf5_dir, f)) and os.path.splitext(f)[1] == '.hdf5']

    hfd5_files = natsorted(hfd5_files)
    pool = Pool()
    _res = pool.map(partial(cell_boundaries_px_helper, cfg), hfd5_files)
    pool.close()
    pool.join()

    # get the results in two lists
    keys = [d[0] for d in _res]
    boundaries = [d[1] for d in _res]

    # flatten the lists
    keys = list(chain.from_iterable(keys))
    boundaries = list(chain.from_iterable(boundaries))

    out = pd.DataFrame({'cell_key': keys,
                        'cell_label': np.arange(1, len(boundaries) + 1).astype(np.int),
                        'cell_boundaries': boundaries})
    # out.to_csv('cell_boundaries.csv', index=False)
    return out


def cell_boundaries_px_helper(cfg, hdf5_file):
    logger.info('rotating all cell boundaries in %s by %d degrees' % (hdf5_file, cfg['rotation'][0]))
    boundaries = []
    keys = []
    fov_cells, fov_cell_keys = get_boundaries(hdf5_file, cfg)
    for fov_cell, fov_cell_key in zip(fov_cells, fov_cell_keys):
        fov_cell = rotate_data(fov_cell, cfg).astype(np.int32)
        boundaries.append(fov_cell.tolist())
        keys.append(fov_cell_key)
    return keys, boundaries



if __name__=="__main__":
    boundaries = cell_boundaries()
    target_dir = os.path.join(config.PROJECT_DIR, 'dashboard', 'data')
    write_tsv(boundaries, target_dir)
    print('Done!')

    # metadata_csv = config.DEFAULT['cell_metadata']
    # hdf5_dir = os.path.join(config.ROOT_DIR, 'cell_boundaries')
    # # cellMetadata = pd.read_csv(metadata_csv).rename(columns={'Unnamed: 0': 'uid'})
    #
    # zIndex = 3
    # fov = 0
    #
    # hfd5_files = [f for f in listdir(hdf5_dir) if isfile(join(hdf5_dir, f)) and os.path.splitext(f)[1] == '.hdf5']
    # out = get_boundaries(fov, zIndex)
    # print(out)
