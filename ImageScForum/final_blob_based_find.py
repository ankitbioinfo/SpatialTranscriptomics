import matplotlib.pyplot as plt
import matplotlib
import xarray as xr
import os
from slicedimage import ImageFormat
from showit import image
import numpy as np
import napari
import seaborn as sns
import starfish.data

from starfish.morphology import Binarize, Filter, Merge, Segment
from starfish import data,Experiment, FieldOfView, display, IntensityTable
from starfish.experiment.builder import write_experiment_json
from starfish.image import Filter as imFilter
from starfish.image import ApplyTransform, LearnTransform
from starfish.types import Axes, Features, Levels, TraceBuildingStrategies
from starfish.util.plot import diagnose_registration,imshow_plane, intensity_histogram
from starfish.spots import DecodeSpots,FindSpots,DetectPixels,DecodeSpots,AssignTargets
from starfish.core.spots.DecodeSpots.trace_builders import build_spot_traces_exact_match,\
build_traces_sequential, build_traces_nearest_neighbors

def segmentCell(imgs,nuclei):
    # set parameters
    dapi_thresh = .18  # global threshold value for nuclei images
    stain_thresh = .22  # global threshold value for primary images
    min_dist = 57  # minimum distance (pixels) between nuclei distance transformed peaks
    min_allowed_size = 10  # minimum size (in pixels) of nuclei
    max_allowed_size = 10000  # maximum size (in pixels) of nuclei

    mp = imgs.reduce({Axes.CH, Axes.ZPLANE}, func="max")
    stain = mp.reduce(
        {Axes.ROUND},
        func="mean",
        level_method=Levels.SCALE_BY_IMAGE)

    nuclei_mp_scaled = nuclei.reduce(
        {Axes.ROUND, Axes.CH, Axes.ZPLANE},
        func="max",
        level_method=Levels.SCALE_BY_IMAGE)

    f = plt.figure(figsize=(12,5))
    ax1 = f.add_subplot(121)
    nuclei_numpy = nuclei_mp_scaled._squeezed_numpy(Axes.ROUND, Axes.CH, Axes.ZPLANE)
    image(nuclei_numpy, ax=ax1, size=20, bar=True)
    plt.title('Nuclei')

    ax2 = f.add_subplot(122)
    image(
        stain._squeezed_numpy(Axes.ROUND, Axes.CH, Axes.ZPLANE),
        ax=ax2, size=20, bar=True)
    plt.title('Stain')
    plt.savefig('1Nuclei_staining.png')
    plt.clf()

    binarized_nuclei = Binarize.ThresholdBinarize(dapi_thresh).run(nuclei_mp_scaled)
    labeled_masks = Filter.MinDistanceLabel(min_dist, 1).run(binarized_nuclei)
    watershed_markers = Filter.AreaFilter(min_area=min_allowed_size, max_area=max_allowed_size).run(labeled_masks)

    plt.subplot(121)
    image(
        binarized_nuclei.uncropped_mask(0).squeeze(Axes.ZPLANE.value).values,
        bar=False,
        ax=plt.gca(),
    )
    plt.title('Nuclei Thresholded')

    plt.subplot(122)
    image(
        watershed_markers.to_label_image().xarray.squeeze(Axes.ZPLANE.value).values,
        size=20,
        cmap=plt.cm.nipy_spectral,
        ax=plt.gca(),
    )
    plt.title('Found: {} cells'.format(len(watershed_markers)))
    plt.savefig('2nuclei_threshold_watershed_marker.png')
    plt.clf()

    thresholded_stain = Binarize.ThresholdBinarize(stain_thresh).run(stain)
    markers_and_stain = Merge.SimpleMerge().run([thresholded_stain, watershed_markers])
    watershed_mask = Filter.Reduce(
        "logical_or",
        lambda shape: np.zeros(shape=shape, dtype=np.bool)
    ).run(markers_and_stain)

    image(
        watershed_mask.to_label_image().xarray.squeeze(Axes.ZPLANE.value).values,
        bar=False,
        ax=plt.gca(),
    )
    plt.title('Watershed Mask')
    plt.savefig('3watershed_mask.png')
    plt.clf()

    segmenter = Segment.WatershedSegment(connectivity=np.ones((1, 3, 3), dtype=np.bool))

    # masks is BinaryMaskCollection for downstream steps
    masks = segmenter.run(
        stain,
        watershed_markers,
        watershed_mask,
    )

    # display result
    image(
        masks.to_label_image().xarray.squeeze(Axes.ZPLANE.value).values,
        size=20,
        cmap=plt.cm.nipy_spectral,
        ax=plt.gca(),
    )
    plt.title('Segmented Cells')
    plt.savefig('4segmentedCells.png')
    plt.clf()
    return segmenter,masks


def plot_intensity_histograms(stack: starfish.ImageStack, r: int,savename: str):
    fig = plt.figure(dpi=150)
    ax1 = fig.add_subplot(231, title='ch: 0')
    ax2 = fig.add_subplot(232, title='ch: 1', sharex=ax1, sharey=ax1)
    ax3 = fig.add_subplot(233, title='ch: 2', sharex=ax1, sharey=ax1)
    ax4 = fig.add_subplot(234, title='ch: 3', sharex=ax1, sharey=ax1)
    ax5 = fig.add_subplot(235, title='ch: 4', sharex=ax1, sharey=ax1)
    ax6 = fig.add_subplot(236, title='ch: 5', sharex=ax1, sharey=ax1)
    intensity_histogram(stack, sel={Axes.ROUND: 0, Axes.CH: 0}, log=True, bins=50, ax=ax1)
    intensity_histogram(stack, sel={Axes.ROUND: 0, Axes.CH: 1}, log=True, bins=50, ax=ax2)
    intensity_histogram(stack, sel={Axes.ROUND: 0, Axes.CH: 2}, log=True, bins=50, ax=ax3)
    intensity_histogram(stack, sel={Axes.ROUND: 1, Axes.CH: 0}, log=True, bins=50, ax=ax4)
    intensity_histogram(stack, sel={Axes.ROUND: 1, Axes.CH: 1}, log=True, bins=50, ax=ax5)
    intensity_histogram(stack, sel={Axes.ROUND: 1, Axes.CH: 2}, log=True, bins=50, ax=ax6)
    fig.tight_layout()
    fig.savefig('Intensity_histogram_'+savename+'.png')


outputdir='Andy_output/'
e = Experiment.from_json(os.path.join(outputdir, "experiment.json"))
RNA_and_nuclei = e.fov().get_image('primary')
nuclei_ref = e.fov().get_image('dots')
nuclei = e.fov().get_image('nuclei')

print('RNA and nuclei channel',RNA_and_nuclei)
print('nuclei reference',nuclei)

imgs=RNA_and_nuclei.sel({Axes.ROUND: (0, 2), Axes.CH: (0,2), Axes.ZPLANE: 0})
all_nuclei=RNA_and_nuclei.sel({Axes.ROUND: (0, 2), Axes.CH: 3, Axes.ZPLANE: 0})
print('nuclei from every rounds',all_nuclei)
print('RNA imgs',imgs)

print(len(e.codebook))


f,(ax1,ax2)=plt.subplots(ncols=2)
matplotlib.rcParams["figure.dpi"] = 250
diagnose_registration(imgs.reduce({Axes.CH}, func="max"), {Axes.ROUND:0}, {Axes.ROUND:1},ax=ax1,title='pre-registered', vmin=0, vmax=255)


from starfish.image import LearnTransform

learn_translation = LearnTransform.Translation(reference_stack=nuclei_ref, axes=Axes.ROUND, upsampling=1000)
transforms_list = learn_translation.run(all_nuclei)


transforms_list.to_json('transforms_list.json')
print(transforms_list)

from starfish.image import ApplyTransform

warp = ApplyTransform.Warp()
#registered_nuclei = warp.run(nucleiRounds, transforms_list=transforms_list, in_place=False)
registered_imgs = warp.run(imgs, transforms_list=transforms_list, in_place=False)


diagnose_registration(registered_imgs.reduce({Axes.CH}, func="max"), {Axes.ROUND:0}, {Axes.ROUND:1},ax=ax2,title='registered')

f.tight_layout()
f.savefig('Registration.png')

#plot_intensity_histograms(stack=imgs, r=0, savename='ori')
cptz_2= imFilter.ClipPercentileToZero(p_min=80, p_max=99.999, level_method=Levels.SCALE_BY_CHUNK)
clipped_both_scaled = cptz_2.run(registered_imgs, in_place=False)
#plot_intensity_histograms(stack=clipped_both_scaled, r=0, savename='ClipPercentileToZero_values_below_80p_and_above_99_999p_and_scale')


def find_spots(imgs):
    p = FindSpots.BlobDetector(
        #min_sigma=(1, 1,1, 1,1,1),
        #max_sigma=(10, 10, 10,10,10,10),
        #num_sigma=30,
        #threshold=0.1,
        #is_volume=False,
        #measurement_type='mean',

        min_sigma=(0.5, 0.5, 0.5),
        max_sigma=(30, 30, 30),
        num_sigma=30,
        threshold=0.1
    )
    # spots from reference
    #spots_from_ref = p.run(image_stack=imgs, reference_image=dots)
    #print('A: Build trace with EXACT_MATCH')
    #print(build_spot_traces_exact_match(spots_from_ref))

    # spots from stack
    spots_from_stack = p.run(image_stack=imgs)
    #print('\nB: Build trace with SEQUENTIAL')
    #print(build_traces_sequential(spots_from_stack))

    '''
    #spots from nearest neighbor with stack
    print('\nC: Build trace with NEAREST_NEIGHBORS')
    print(build_traces_nearest_neighbors(spots_from_stack, search_radius=5))

    #spots from nearest neighbor with reference
    print('\nD: Build trace with NEAREST_NEIGHBORS')
    print(build_traces_nearest_neighbors(spots_from_ref, search_radius=5))
    '''

    #intensities=spots_from_ref
    return spots_from_stack


seg,masks=segmentCell(clipped_both_scaled,nuclei)
spots = find_spots(clipped_both_scaled)


def decode_spots(codebook, spots):
    #decoder = DecodeSpots.MetricDistance(codebook=codebook,max_distance=1,min_intensity=1,metric='euclidean',
    #           norm_order=2,trace_building_strategy=TraceBuildingStrategies.EXACT_MATCH)

    #decoder = DecodeSpots.MetricDistance(codebook=codebook,max_distance=0.2,min_intensity=0.1,
    #           trace_building_strategy=TraceBuildingStrategies.NEAREST_NEIGHBOR)

    #decoder = DecodeSpots.MetricDistance(codebook=codebook,max_distance=0.2,min_intensity=0.1,metric='euclidean',
    #                                     norm_order=2,trace_building_strategy=TraceBuildingStrategies.EXACT_MATCH)

    decoder= DecodeSpots.PerRoundMaxChannel(codebook=codebook,anchor_round=0,search_radius=100,
                                  trace_building_strategy=TraceBuildingStrategies.NEAREST_NEIGHBOR)

    #decoder= DecodeSpots.PerRoundMaxChannel(codebook=codebook,anchor_round=0,search_radius=100,
    #                              trace_building_strategy=TraceBuildingStrategies.EXACT_MATCH)


    #decoded = initial_spot_intensities.loc[spots[Features.PASSES_THRESHOLDS]]
    #decoded_filtered = decoded[decoded.target != 'nan']
    #decoder = DecodeSpots.SimpleLookupDecoder(codebook=codebook)
    return decoder.run(spots=spots)

def make_expression_matrix(masks, decoded):
    #decoded_filtered = decoded["fov_000"].loc[decoded["fov_000"][Features.PASSES_THRESHOLDS]]
    al = AssignTargets.Label()
    labeled = al.run(masks, decoded[decoded.target != 'nan'])
    labeled_filtered = labeled[labeled.cell_id != 'nan']
    cg=labeled_filtered.to_expression_matrix()
    return labeled_filtered,cg


decoded = decode_spots(e.codebook, spots)
genes, counts = np.unique(decoded.loc[decoded[Features.PASSES_THRESHOLDS]][Features.TARGET], return_counts=True)
print('genes and counts', len(genes),len(counts))
for i in range(len(genes)):
    print(genes[i],'\t',counts[i])


labeled_filtered,mat = make_expression_matrix(masks, decoded)
genes=mat.coords['genes'].values
X=mat.coords['x'].values
Y=mat.coords['y'].values
for i in range(len(genes)):
    print(X[i],Y[i],genes[i])

print(mat.dims)
print(mat.values.shape)
print([genes,np.sum(mat.values,axis=0)] )
