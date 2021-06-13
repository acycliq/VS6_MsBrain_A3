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
import src.config as config
from src.utils import splitter_mb, transformation
import logging

logger = logging.getLogger()
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s:%(levelname)s:%(message)s"
)


def write_boundaries_tsv(boundaries, target_dir):
    N = len(boundaries)
    dummy_list = [[] for d in boundaries]
    dummy_prob = [[1.0] for d in boundaries]
    dummy_name = [['Generic'] for d in boundaries]
    df = pd.DataFrame({'Cell_Num': 1 + np.arange(N),
                       'X': np.zeros(N).astype(np.int),
                       'Y': np.zeros(N).astype(np.int),
                       'Genenames': dummy_list,
                       'CellGeneCount': dummy_list,
                       'ClassName': dummy_name,
                       'Prob': dummy_prob
                       })
    df.to_csv(os.path.join(target_dir, 'cellData.tsv'), sep='\t', index=False)
    splitter_mb(df, os.path.join(target_dir, 'cellData'), 99)


def write_celldata_tsv(boundaries, target_dir):
    df = pd.DataFrame({'cell_id': 1+ np.arange(len(boundaries)),
                       'coords': boundaries})
    df.to_csv(os.path.join(target_dir, 'cellBoundaries.tsv'), sep='\t', index=False)
    splitter_mb(df, os.path.join(target_dir, 'cellBoundaries'), 99)


def write_tsv(boundaries, out_dir):
    write_boundaries_tsv(boundaries, out_dir)
    write_celldata_tsv(boundaries, out_dir)


def micron_to_pixel(temp, cfg):
    tx, ty, _, _ = transformation(cfg)
    x = tx(temp[:, 0]).astype(np.int)
    y = ty(temp[:, 1]).astype(np.int)
    x_max = x.max()
    y_max = y.max()
    return list(zip(x, y)), x_max, y_max


def get_boundaries(fileName: str, cfg, zIndex:int=3):
    """
    Helper function that reads a cached hdf5 file and gets the cell boundaries for that specific fov only
    :param fileName:
    :param zIndex:
    :return:
    """
    target_hdf5 = os.path.join(config.ROOT_DIR, 'cell_boundaries', fileName)
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
    hdf5_dir = os.path.join(config.ROOT_DIR, 'cell_boundaries')
    # metadata_csv = config.DEFAULT['cell_metadata']
    # cellMetadata = pd.read_csv(metadata_csv).rename(columns={'Unnamed: 0': 'uid'})

    hfd5_files = [f for f in listdir(hdf5_dir) if isfile(join(hdf5_dir, f)) and os.path.splitext(f)[1] == '.hdf5']
    boundaries = []
    keys = []
    for hdf5_file in natsorted(hfd5_files):
        fov_cells, fov_cell_keys = get_boundaries(hdf5_file, cfg)
        for fov_cell, fov_cell_key in zip(fov_cells, fov_cell_keys):
            boundaries.append(fov_cell)
            keys.append(fov_cell_key)
    out = pd.DataFrame({'cell_key': keys,
                        'cell_label': np.arange(1, len(boundaries) + 1).astype(np.int),
                        'cell_boundaries': boundaries})
    out.to_csv('cell_boundaries.csv', index=False)
    return out


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
