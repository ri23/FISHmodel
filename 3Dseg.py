"""Segment 3D tissue without cell walls."""

import os
import argparse

import numpy as np

import scipy.ndimage
import scipy.misc
from scipy.ndimage.filters import laplace

from skimage.exposure import equalize_hist
from skimage.filters import gaussian_filter
from skimage.measure import label
from skimage.morphology import watershed, remove_small_objects

from jicbioimage.core.io import FileBackend
from jicbioimage.core.image import DataManager
from jicbioimage.core.image import SegmentedImage

from jicbioimage.transform import (
    max_intensity_projection
)

from jicbioimage.illustrate import AnnotatedImage

HERE = os.path.dirname(__file__)
UNPACK = os.path.join(HERE, '..', 'data', 'unpack')
OUTPUT = os.path.join(HERE, '..', 'output')#'/group-share','ietswaar','test','output')#HERE, '..', 'output') RI edit 1

if not os.path.isdir(OUTPUT):
    os.mkdir(OUTPUT)

DEBUG = False

def collection_from_filename(stack_filename):
    file_backend = FileBackend(UNPACK)
    data_manager = DataManager(file_backend)
    microscopy_collection = data_manager.load(stack_filename)

    return microscopy_collection

def save_sample(filename, stack, sample_z=25):

    full_path = os.path.join(OUTPUT, filename)
    if DEBUG:
        scipy.misc.imsave(full_path, stack[:,:,sample_z])

def save_stack(stack, stack_name='stack'):

    if not DEBUG:
        return

    stack_dir = os.path.join(OUTPUT, stack_name + '.stack')

    if not os.path.isdir(stack_dir):
        os.mkdir(stack_dir)

    xdim, ydim, zdim = stack.shape

    for z in range(zdim):
        filename = 'z{}.png'.format(z)
        full_name = os.path.join(stack_dir, filename)
        scipy.misc.imsave(full_name, stack[:,:,z])

def blank_layers(input_array, n_layers=2, blank=1):
    """Return a copy of the input array with the top and bottom
    n_layers set to a particular value."""

    _, _, zdim = input_array.shape

    start_z = n_layers
    stop_z = zdim - n_layers

    blanked = input_array.copy()

    blanked[:,:,0:start_z] = blank
    blanked[:,:,stop_z:] = blank

    return blanked

def find_seeds(zstack):
    """Return array containing segmentation seeds."""

    smooth_sigma = 10
    seed_threshold = 0.13
    min_size = 40000#10000 RI edit 5

    xdim, ydim, zdim = zstack.shape

    save_sample('start.png', zstack)

    smoothed = gaussian_filter(zstack, sigma=smooth_sigma)
    save_sample('smoothed.png', smoothed)

    edges = laplace(smoothed)
    edges = edges + np.min(edges)
    save_sample('laplace.png', edges)

    equalised = equalize_hist(edges)
    save_sample('equalised.png', equalised)

    blanked = blank_layers(equalised)
    thresholded = blanked < seed_threshold

    save_sample('thresholded.png', thresholded)
    save_stack(thresholded, 'thresh')

    connected = label(thresholded)
    save_sample('connected.png', connected)
    save_stack(connected, 'connected')

    #rids = np.unique(connected)
    #print [len(np.where(connected==rid)[0]) for rid in rids[1:]]

    filtered_connected = remove_small_objects(connected, min_size=min_size)
    save_stack(filtered_connected, 'filtered_connected')

    return filtered_connected

def segment_from_seeds(zstack, seeds, watershed_cutoff):

    smooth_sigma =5 #15 RI edit 4
    size_threshold = 10000

    smoothed2 = scipy.ndimage.filters.gaussian_filter(zstack, 
                                                      sigma=smooth_sigma)
    save_sample('smoothed2.png', smoothed2)

    inverted = np.max(smoothed2) - smoothed2
    save_sample('inverted.png', inverted)

    # Now normalised
    equalised2 = equalize_hist(inverted)
    save_sample('equalised2.png', equalised2)
    save_stack(equalised2, 'equalised')

    mask = equalised2 < watershed_cutoff
    save_sample('mask.png', mask)

    segmented = watershed(equalised2, seeds, mask=mask)
    save_sample('segmented.png', segmented)
    save_stack(segmented, 'segmented')

    # region_ids = np.unique(segmented)
    # sizes = [len(np.where(segmented == rid)[0]) for rid in region_ids]

    nosmall = remove_small_objects(segmented, min_size=size_threshold)
    save_stack(nosmall, 'nosmall')

    reseg = watershed(equalised2, nosmall, mask=mask)
    save_stack(reseg, 'reseg')

    return reseg

def uint8ify(input_array):

    max_val = float(np.max(input_array))
    min_val = float(np.min(input_array))

    val_range = max_val - min_val

    return 255 * ((input_array.astype(np.float) - min_val) / val_range)

def generate_annotated_image(collection, cell_level_threshold):

    zstack = collection.zstack_array(s=0, c=2)
    probe_stack = collection.zstack_array(s=0, c=0)

    max_intensity_projection(probe_stack)

    seeds = find_seeds(zstack)

    #probe_stack2 = collection.zstack_array(s=0, c=1) #RI edit 2
    zstack = zstack + probe_stack #+ probe_stack2#RI edit 3

    segmentation = segment_from_seeds(zstack, seeds, cell_level_threshold)

    projection = max_intensity_projection(zstack)
    projection_as_uint8 = uint8ify(projection)
    annotated_projection = AnnotatedImage.from_grayscale(projection_as_uint8)

    
    rids = np.unique(segmentation)

    for rid in rids[1:]:
        x, y, z = map(np.mean, np.where(segmentation == rid))
        size = len(np.where(segmentation == rid)[0])
    
        annotated_projection.text_at(str(size), y-10, x)

    annotation_filename = 'annotated_image.png'
    with open(annotation_filename, 'wb') as f:
        f.write(annotated_projection.png())

def main():

    global DEBUG

    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('image_filename', help="Image filename")
    parser.add_argument('--cell-level-threshold', 
                        type=float,
                        default=0.3,
                        help="Threshold (in range 0 < t < 1) defining cell")
    parser.add_argument('--verbose',
                        type=bool,
                        default=False,
                        help="Whether processing stages should be output")
    args = parser.parse_args()

    DEBUG = args.verbose

    collection = collection_from_filename(args.image_filename)

    generate_annotated_image(collection, args.cell_level_threshold)

if __name__ == "__main__":
    main()
