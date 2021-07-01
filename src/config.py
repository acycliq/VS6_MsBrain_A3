import os
from pathlib import Path
from collections import defaultdict

PROJECT_DIR = Path(os.path.dirname(os.path.realpath(__file__))).parent.absolute()
print(PROJECT_DIR)

DROPBOX_URL = os.path.join('Z:\\', 'MERFISH_F_E')
def get_config(slice_id, region_id):
    out = {
        'cell_by_gene': os.path.join(DROPBOX_URL, slice_id, region_id, 'cell_by_gene.csv'),
        'cell_metadata': os.path.join(DROPBOX_URL, slice_id, region_id, 'cell_metadata.csv'),
        'detected_transcripts': os.path.join(DROPBOX_URL, slice_id, region_id, 'detected_transcripts.csv'),
        # 'manifest': os.path.join(DROPBOX_URL, slice_id, region_id, 'images', 'manifest.json'),
        'micron_to_mosaic_pixel_transform': os.path.join(DROPBOX_URL, slice_id, region_id, 'images', 'micron_to_mosaic_pixel_transform.csv'),
        'dapi_tif': os.path.join(DROPBOX_URL, slice_id, region_id, 'images', 'mosaic_DAPI_z3.tif'),
        'cell_boundaries_dir': os.path.join(DROPBOX_URL, slice_id, region_id, 'cell_boundaries'),
        'clip_poly': _clip_poly[slice_id][region_id],
        'rotation': _rotation[slice_id][region_id],
    }
    return out


# only spots within this poly will be plotted. Use empty array if you do not want any clipping to be applied
_clip_poly = defaultdict(dict)
_clip_poly['MsBrain_Eg1_VS6_JH_V6_05-02-2021']['region_0'] = [(23000, 30684), (24323, 51017), (72822, 39953), (64364, 23533), (43641, 28538)]
_clip_poly['MsBrain_Eg1_VS6_JH_V6_05-02-2021']['region_1'] = [(23080, 29550), (24566, 47472), (70552, 37593), (64694, 24042), (43800, 27451)]

_clip_poly['MsBrain_Eg2_VS6_V11_JH_05-02-2021']['region_0'] = [(25740, 24217), (12794, 33203), (47291, 60999), (55363, 45387)]
_clip_poly['MsBrain_Eg2_VS6_V11_JH_05-02-2021']['region_1'] = [(39595, 33354), (48856, 18552), (18254, -1342), (10335, 9731)]

_clip_poly['MsBrain_Eg3_VS6_JH_V6_05-01-2021']['region_0'] = [(25255, 29712), (27483, 35489), (27648, 40028), (27731, 41596), (26245, 44072), (24430, 47043), (23192, 47951), (23934, 55094), (32353, 66356), (41431, 44898), (38708, 23109), (30950, 26410), (25750, 29712)]

# just guessing the clipping poly...
_clip_poly['MsBrain_Eg3_VS6_JH_V6_05-01-2021']['region_1'] = [(25255, 29712), (27483, 35489), (27648, 40028), (27731, 41596), (26245, 44072), (24430, 47043), (23192, 47951), (23934, 55094), (32353, 66356), (41431, 44898), (38708, 23109), (30950, 26410), (25750, 29712)]

_clip_poly['MsBrain_EG4_VS6library_V6_LH_04-14-21']['region_0'] = [(35367, 40796), (35485, 47248), (54566, 47563), (54566, 41072)]

# just guessing the clipping poly...
_clip_poly['MsBrain_EG4_VS6library_V6_LH_04-14-21']['region_1'] = [(35367, 40796), (35485, 47248), (54566, 47563), (54566, 41072)]

_clip_poly['MsBrain_Eg5_VS6_JH_V6_05-16-2021']['region_0'] = [(25493, 20300), (25650, 31158), (54841, 32496), (54943, 18648)]
_clip_poly['MsBrain_Eg5_VS6_JH_V6_05-16-2021']['region_1'] = [(24674, 16221), (25968, 26958), (53383, 25511), (54221, 15535)]


_rotation = defaultdict(dict)
_rotation["MsBrain_Eg1_VS6_JH_V6_05-02-2021"]["region_0"] = 9,
_rotation["MsBrain_Eg1_VS6_JH_V6_05-02-2021"]["region_1"] = 11,
_rotation["MsBrain_Eg2_VS6_V11_JH_05-02-2021"]["region_0"] = -39,
_rotation["MsBrain_Eg2_VS6_V11_JH_05-02-2021"]["region_1"] = -42,
_rotation["MsBrain_Eg3_VS6_JH_V6_05-01-2021"]["region_0"] = 75,
_rotation["MsBrain_Eg3_VS6_JH_V6_05-01-2021"]["region_1"] = 79,
_rotation['MsBrain_EG4_VS6library_V6_LH_04-14-21']['region_0'] = 0,
_rotation['MsBrain_EG4_VS6library_V6_LH_04-14-21']['region_1'] = 0,
_rotation['MsBrain_Eg5_VS6_JH_V6_05-16-2021']['region_0'] = 0,
_rotation['MsBrain_Eg5_VS6_JH_V6_05-16-2021']['region_1'] = 0,
