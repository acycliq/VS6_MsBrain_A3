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
_clip_poly['MsBrain_Eg1_VS6_JH_V6_05-02-2021']['region_0'] = [(25009,29875), (25594,32860), (27268,37773), (28476,42200), (28319,45239), (48535,40838), (72239,39897), (71366,38387), (68110,34895), (67266,32498), (63527,25493), (52594,26965), (45625,27624), (33621,28185)]
_clip_poly['MsBrain_Eg1_VS6_JH_V6_05-02-2021']['region_1'] = [(25106,29454), (26388,32536), (26659,38544), (28245,41168), (27516,44897), (48331,39522), (69056,36621), (65948,33412), (66650,32191), (65446,28944), (65549,22916), (45505,28184)]

_clip_poly['MsBrain_Eg2_VS6_V11_JH_05-02-2021']['region_0'] = [(28094,23896), (25746,25410), (24495,27027), (20162,29651), (17683,30951), (15614,31536), (33310,44154), (50425,59685), (51163,54493), (52558,47763), (55733,45640), (58030,44952), (47527,38343), (43377,38534), (40182,36575), (38295,34134), (38376,30361)]
_clip_poly['MsBrain_Eg2_VS6_V11_JH_05-02-2021']['region_1'] = [(22166,21044), (40780,31079), (43060,25614), (44783,20500), (48653,17247), (51944,17096), (39781,8273), (36062,8532), (33925,7331), (32613,6355)]

_clip_poly['MsBrain_Eg3_VS6_JH_V6_05-01-2021']['region_0'] = [(19359,45462), (21169,57544), (23716,57658), (24839,58191), (27401,60449), (30586,63612), (32761,66082), (43803,45026), (41638,22295), (30945,26491), (29633,26460), (28078,25179), (26332,23426), (22038,30291), (23938,30075), (25134,29845), (26418,31767), (26731,33919), (27215,37128), (27286,41473), (26035,43890), (24784,45855), (23503,46923), (21318,46846)]
_clip_poly['MsBrain_Eg3_VS6_JH_V6_05-01-2021']['region_1'] = [(20118,47441), (24203,54149), (26177,56640), (30439,61406), (40574,41546), (39264,20452), (37252,20928), (35823,20936), (27965,25666), (25275,25646), (27374,39593)]

_clip_poly['MsBrain_EG4_VS6library_V6_LH_04-14-21']['region_0'] = [(38795,34051), (48315,36706), (53254,35519), (57565,33713), (60644,28831), (63749,26289), (39603,18535), (16894,23371), (19447,28288), (21805,31870), (24055,34746), (27814,35267), (32734,34724)]
_clip_poly['MsBrain_EG4_VS6library_V6_LH_04-14-21']['region_1'] = [(39155,31287), (40553,31074), (42232,31653), (45079,32768), (48324,33984), (50624,34889), (53553,34110), (56277,32898), (58316,31663), (59659,30535), (64268,25448), (42512,16254), (19976,19858), (20980,21404), (23078,25355), (26402,29978), (28978,31853), (31479,31942), (34826,31488), (37966,31032)]

_clip_poly['MsBrain_Eg5_VS6_JH_V6_05-16-2021']['region_0'] = [(39035,30314), (45997,32143), (48813,33976), (51279,33900), (55589,33350), (59792,31342), (63320,25053), (40458,15638), (17272,20468), (17414,23563), (16671,25882), (16976,27139), (17830,28421), (20984,30671), (24045,31357), (27619,31120), (31897,31233), (36821,30905)]
_clip_poly['MsBrain_Eg5_VS6_JH_V6_05-16-2021']['region_1'] = [(38880,25636), (42937,25977), (46003,27877), (48962,29798), (52568,29600), (54399,28868), (57815,26794), (63062,21917), (63992,17897), (39511,10724), (15703,14987), (16212,19046), (17859,22727), (20355,24649), (22694,25969), (27225,27157), (30450,27080), (33618,26570)]


_rotation = defaultdict(dict)
_rotation["MsBrain_Eg1_VS6_JH_V6_05-02-2021"]["region_0"] = 9,
_rotation["MsBrain_Eg1_VS6_JH_V6_05-02-2021"]["region_1"] = 11,
_rotation["MsBrain_Eg2_VS6_V11_JH_05-02-2021"]["region_0"] = -39,
_rotation["MsBrain_Eg2_VS6_V11_JH_05-02-2021"]["region_1"] = -42,
_rotation["MsBrain_Eg3_VS6_JH_V6_05-01-2021"]["region_0"] = 75,
_rotation["MsBrain_Eg3_VS6_JH_V6_05-01-2021"]["region_1"] = 79,
_rotation['MsBrain_EG4_VS6library_V6_LH_04-14-21']['region_0'] = 180,
_rotation['MsBrain_EG4_VS6library_V6_LH_04-14-21']['region_1'] = 180,
_rotation['MsBrain_Eg5_VS6_JH_V6_05-16-2021']['region_0'] = 180,
_rotation['MsBrain_Eg5_VS6_JH_V6_05-16-2021']['region_1'] = 180,
