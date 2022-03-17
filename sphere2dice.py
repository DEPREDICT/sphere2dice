# This code flattens 3D surface data in to a selected number of patches.
#
# M.G.Poirot November 2021. Refactored in February 2022.

import nibabel.freesurfer.io as fsio
import matplotlib.pyplot as plt
import numpy as np
import nibabel as nib
from scipy.spatial.transform import Rotation as R
from scipy.interpolate import griddata
from os.path import join, basename, isfile
from os import makedirs
from tqdm import tqdm
from copy import deepcopy
import sys
import argparse
from os.path import splitext

from collections import OrderedDict

# MutableSet class was moved after Python 3.3, so for backward compatibility:
try:
    from collections.abc import MutableSet
except ModuleNotFoundError or ImportError:
    from collections import MutableSet


class OrderedSet(MutableSet):
    '''
        Ordered sets are not a Python Builtin. 
        We need it to keep track of our rotations 
        whilst keeping them in order.

        Do not be bothered by this bulky code too much
        it is as simple as it sounds.
    '''

    def __init__(self, iterable=None):
        self.end = end = []
        end += [None, end, end]  # sentinel node for doubly linked list
        self.map = {}  # key --> [key, prev, next]
        if iterable is not None:
            self |= iterable

    def __len__(self):
        return len(self.map)

    def __contains__(self, key):
        return key in self.map

    def add(self, key):
        if key not in self.map:
            end = self.end
            curr = end[1]
            curr[2] = end[1] = self.map[key] = [key, curr, end]

    def discard(self, key):
        if key in self.map:
            key, prev, next = self.map.pop(key)
            prev[2] = next
            next[1] = prev

    def __iter__(self):
        end = self.end
        curr = end[2]
        while curr is not end:
            yield curr[0]
            curr = curr[2]

    def __reversed__(self):
        end = self.end
        curr = end[1]
        while curr is not end:
            yield curr[0]
            curr = curr[1]

    def pop(self, last=True):
        if not self:
            raise KeyError('set is empty')
        key = self.end[1][0] if last else self.end[2][0]
        self.discard(key)
        return key

    def __repr__(self):
        if not self:
            return '%s()' % (self.__class__.__name__,)
        return '%s(%r)' % (self.__class__.__name__, list(self))

    def __eq__(self, other):
        if isinstance(other, OrderedSet):
            return len(self) == len(other) and list(self) == list(other)
        return set(self) == set(other)


def _rotate_sphere(xz_turn: tuple, _is_right: bool, _sphere: np.ndarray):
    """
    Return rotated sphere coordinates as rotated in Z-X order,
    taking into account flipping Right Hemispheres over X-axis to match left hemi.
    :param xz_turn: A tuple of length 2 with X and Z rotation respectively
    :param _is_right: A boolean
    :param _sphere: A numpy array containing sphere coordinates
    :return sphere_r: A numpy array containing rotated sphere coordinates
    """
    # Extract and copy variables
    x_turn, z_turn = xz_turn
    rotated_sphere = deepcopy(_sphere)
    # Flip over X-axis if if is a right hemisphere
    if _is_right:
        # Flip axis for right to match left
        rotated_sphere[:, 0] *= -1
        # We found out that we need to correct the right hemisphere for a slight wobble
        _correction_turn = (-4.11, 12.89, -10.56)  # USED CORRECTION FROM GRIDSEARCH
        _correction = R.from_rotvec(np.radians(np.array(_correction_turn)))
        rotated_sphere = _correction.apply(rotated_sphere)
    # Rotate over Z-axis if the turn over Z is not 0
    if z_turn:
        _rotation = R.from_rotvec(np.radians(np.array((0, 0, z_turn))))
        rotated_sphere = _rotation.apply(rotated_sphere)
    # Rotate over X-axis if the turn over X is not zero
    if x_turn:
        _rotation = R.from_rotvec(np.radians(np.array((x_turn, 0, 0))))
        rotated_sphere = _rotation.apply(rotated_sphere)
    return rotated_sphere


def _process_fig(create_patch, patch, scalars, _is_right, xz_rotation,
                 patch_labels=None):
    """
    Plotting this figure is just for verbosity/debugging purposes.
    It creates a 3x5 subplot of 5 figures displaying the generation process.
    """

    def _norm_axes(ax, _is_right: False):
        """ Applies default formatting for 3D objects
        :param ax: matplotlib.pyplot.axis object
        :param _is_right: boolean used for more accurate axis labeling
        """
        # Note that we flipped the coordinates
        flipped = ''
        if _is_right:
            flipped = '-'

        ax.set_xlabel(flipped + 'x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        ax.set_xlim([-1, 1])
        ax.set_ylim([-1, 1])
        ax.set_zlim([-1, 1])
        ax.set_aspect('auto')

    # Retrieve data stored in wrapper
    spheres_r, spheres_f, s = create_patch()
    # Get neutral sphere
    x, y, z = np.split(spheres_r[(False, (0, 0))], 3, axis=1)
    # Get rotated sphere
    d, x_d, z_d = spheres_f[(_is_right, xz_rotation)]
    _n_patches = len(set(rotations for _, rotations in spheres_r))

    # Keep track of which data point corresponds with which patch
    patch_bool = np.logical_and(np.abs(x_d) <= s, np.abs(z_d) <= s)
    if patch_labels is None:
        # On first loop, convert current label to integer
        patch_labels = patch_bool.astype(int)
    else:
        # On subsequent loops, add current label to integer list.
        patch_labels[patch_bool] = np.max(patch_labels) + 1

    # Subplot 1: Original spherical data
    fig = plt.figure(figsize=(12, 8))
    ax1 = fig.add_subplot(2, 3, 1, projection='3d')
    ax1.scatter(xs=x, ys=y, zs=z, c=scalars, s=1, cmap='jet')
    _norm_axes(ax1, _is_right)
    ax1.set_title('Scalar sphere')

    # Subplot 2: Distance to sampling point
    ax2 = fig.add_subplot(2, 3, 2, projection='3d')
    sc1 = ax2.scatter(xs=x, ys=y, zs=z, c=d, s=1, cmap='inferno_r')
    _norm_axes(ax2, _is_right)
    ax2.set_title('Distance to patch center')
    fig.colorbar(sc1)

    # Subplot 3: 2D projected distance to sampling point
    ax3 = fig.add_subplot(2, 3, 3)
    ax3.scatter(x_d, z_d, c=scalars, s=1, cmap='jet')
    ax3.plot([-s, -s, s, s, -s], [-s, s, s, -s, -s], c='k')
    ax3.set_title('Distance to patch center and patch border')

    # Subplot 4: Selected patch on sphere
    ax4 = fig.add_subplot(2, 3, 4, projection='3d')
    ax4.scatter(xs=x, ys=y, zs=z, s=1, c=patch_labels,
                cmap='gist_stern', vmin=0, vmax=_n_patches)
    _norm_axes(ax4, _is_right)
    ax4.set_title('Patch domain')

    # Subplot 5: Content of captured patch
    ax5 = fig.add_subplot(2, 3, 5)
    ax5.imshow(patch, cmap='jet')
    ax5.set_title('Final {}x{} patch'.format(len(patch), len(patch)))

    fig.suptitle('Direction ' + str(xz_rotation) + ' from "{}"')

    # Return the figure handle, and patch labels assigned to each coordinate
    return fig, patch_labels


def get_rotated_spheres_dict(_sphere_file: str, _n_patches=6):
    """
    Get the X and Z rotations corresponding to the desired number of patches.
    :param _n_patches:
    :param _sphere_file:
    :param n_patches: number of patches as integer
    :return: list of tuples of length two, containing X and Z rotation.
    """
    # Load coordinates of a sphere
    sphere, _ = fsio.read_geometry(_sphere_file)

    # Normalize sphere
    sphere = np.subtract(sphere, np.mean(sphere, axis=0))
    dists = np.sqrt(np.sum(sphere ** 2, axis=1))[:, np.newaxis]
    sphere = np.multiply(sphere, np.repeat(1 / dists, 3, axis=1))

    # Define the rotations required for the n_patches options
    # You can add alternate rotation options here, if desired.
    two = [(0, 0), (0, 90)]
    six = two + [(0, 180), (0, -90),
                 (90, 0), (-90, 0)]
    eighteen = six + [(0, 45), (0, 135), (0, 225), (0, 315),
                      (45, 0), (45, 90), (45, 180), (45, 270),
                      (-45, 0), (-45, 90), (-45, 180), (-45, 270)]
    try:
        rotations = {2: two, 6: six, 18: eighteen}[_n_patches]
    except KeyError:
        raise ValueError(f'Ill defined number of patches. Expected 2, 6 or 18, received:{_n_patches}')
    spheres_r = OrderedDict()
    spheres_f = OrderedDict()
    for is_right in range(2):
        for rotation in rotations:
            # Rotate sphere
            sphere_r = _rotate_sphere(rotation, bool(is_right), sphere)
            # Flatten the sphere from pole
            x_r, y_r, z_r = [sphere_r[:, ax] for ax in range(3)]
            # Compute distance over the sphere
            d = np.arccos(-y_r)
            # Extend unit vectors to spherical distance and decompose
            v = np.stack([x_r, z_r])
            x_d, z_d = v / ((v ** 2).sum(axis=0) ** .5) * d
            sphere_f = (d, x_d, z_d)
            # Assign
            spheres_r[bool(is_right), rotation] = sphere_r
            spheres_f[bool(is_right), rotation] = sphere_f
    return spheres_r, spheres_f


def stack_fig(_patch_stack: np.ndarray, source_name='[no file name provided]',
              rotations=None):
    """
    Create a figure with six sides of a dice:
    :param rotations:
    :param _patch_stack: numpy array of size (patch_n, resolution, resolution)
    :param source_name: optional file name used in figure title
    :return: figure handle used for showing or saving
    """
    if rotations is None:
        rotations = 6 * ['']
    else:
        rotations = list(rotations)
    r = _patch_stack.shape[-1]
    stack_figure, ax = plt.subplots()

    # plot six images, each for each side of the dice
    img = np.empty([r * 3, r * 4]) * np.nan
    img[r:r * 2, r * 1:r * 2] = np.flipud(_patch_stack[0])  # center
    img[r:r * 2, r * 0:r * 1] = np.flipud(_patch_stack[1])  # left
    img[r:r * 2, r * 3:r * 4] = np.flipud(_patch_stack[2])  # far
    img[r:r * 2, r * 2:r * 3] = np.flipud(_patch_stack[3])  # right
    img[r * 0:r * 1, r:r * 2] = np.flipud(_patch_stack[4])  # up
    img[r * 2:r * 3, r:r * 2] = np.flipud(_patch_stack[5])  # down
    ax.imshow(img, cmap='jet', vmin=0, vmax=5)

    # place text
    ax.text(r * 1, r, '[0] center ' + str(rotations[0]), va='top')  # center
    ax.text(r * 2, r, '[1] right ' + str(rotations[1]), va='top')  # right
    ax.text(r * 3, r, '[2] far ' + str(rotations[2]), va='top')  # far
    ax.text(r * 0, r, '[3] left ' + str(rotations[3]), va='top')  # left
    ax.text(r * 1, 0, '[4] up ' + str(rotations[4]), va='top')  # up
    ax.text(r * 1, r * 2, '[5] down ' + str(rotations[5]), va='top')  # down

    # formatting
    ax.axis('off')
    ax.axis('equal')
    ax.set_title(f'UV map of 6/{len(_patch_stack)}patches for\n{basename(source_name)}')
    return stack_figure


def make_patch_wrapper(_spheres: tuple, _patch_size=.5 ** .5,
                       _resolution=128, method='linear'):
    spheres_r, spheres_f = _spheres
    # Define interpolation grid
    li = np.linspace(-_patch_size, _patch_size, _resolution)
    ui, vi = np.meshgrid(li, li)

    def _make_patch(scalars=None, is_right=False, xz_rotation=(0, 0)):
        if scalars is None:
            return spheres_r, spheres_f, _patch_size
        d, x_d, z_d = spheres_f[(is_right, xz_rotation)]
        patch = griddata((x_d, z_d), scalars, (ui, vi), method=method)
        patch = np.nan_to_num(patch)
        return patch

    return _make_patch


def make_stack(_make_patch, _scalar_file: str, verbose=False, _is_right=False):
    rotations = OrderedSet(rotations for _, rotations in _make_patch()[0])
    _patch_stack = None
    if splitext(_scalar_file)[-1] == '.gii':
        scalars = nib.load(_scalar_file).agg_data()
    else:
        scalars = fsio.read_morph_data(_scalar_file)
    process_figures = []
    patch_labels = None
    for n, xz_rotation in enumerate(tqdm(rotations)):
        patch = _make_patch(scalars, _is_right, xz_rotation)
        if _patch_stack is None:
            _patch_stack = np.zeros([len(rotations), *patch.shape])
        _patch_stack[n] = patch
        if verbose:
            fig, patch_labels = _process_fig(_make_patch, patch, scalars,
                                             _is_right, xz_rotation,
                                             patch_labels=patch_labels)
            process_figures.append(fig)
    return _patch_stack, process_figures


if __name__ == '__main__':
    # All arguments are optional
    parser = argparse.ArgumentParser(
        description='Turn 3D scalar data from .gii format into a stack of 2D projections in .npy format.')
    parser.add_argument('--input_dirs', type=str, required=False, default=['input'], nargs='*',
                        help='Optional directory to find input files (sphere and scalar files).')
    parser.add_argument('--sphere_file', type=str, required=False, default='fsaverage.sphere',
                        help='generic file with coordinates to which scalars are mapped.')
    parser.add_argument('--scalar_files', type=str, required=False, default=['example-thickness-file.gii'], nargs='*',
                        help='This is where your scalar data goes in. Default is a dummy file.')
    parser.add_argument('--is_rights', type=bool, required=False, default=[False], nargs='*',
                        help='if the hemisphere the data is coming from is the right hemisphere. This will flip and '
                             'perform a slight rotation such that Left and Right hemisphere data line up perfectly.')
    parser.add_argument('--output_dirs', type=str, required=False, default=['output'], nargs='*',
                        help='Optional directory to which output is saved (Like the produced 2D stacks and figures).')
    parser.add_argument('--stack_name', type=str, required=False, default=None,
                        help='filename of the produced 2D stack. Default uses input file name to create one.')
    parser.add_argument('--boxfg_name', type=str, required=False, default=None,
                        help='filename of the graphic of the 2D stack. Default uses input file name to create one.')
    parser.add_argument('--prcfg_name', type=str, required=False, default='angle-{}_overview.jpg',
                        help='filename of an optional set of figures produced when verbose=True. Note that brackets are'
                             ' required to store the appropriate order')
    parser.add_argument('--verbose', required=False, default=False, action='store_true',
                        help='Wether to generate and save figures showing the process of 2D projection.')
    parser.add_argument('--n_patches', type=int, required=False, default=6,
                        help='Number of patches to sample from sphere: options are 2, 6 or 18.')
    parser.add_argument('--resolution', type=int, required=False, default=128,
                        help='resolution of the resulting 2D patch.')
    parser.add_argument('--patch_size', type=float, required=False, default=.5 ** .5,
                        help='rpatch size, make it larger to cover a larger patch of the sphere')
    parser.add_argument('--interpolation', type=str, required=False, default='linear',
                        help='interpolation method for rasterization of coordinates')
    parser.add_argument('--overwrite', type=bool, required=False, default=False,
                        help='boolean to prevent accidental overwriting')
    args, _ = parser.parse_known_args()

    # Construct input
    in_sphere_file = join('input', args.sphere_file)

    # Precompute the rotations and store it in the wrapper. This speeds speeds up the for loop.
    spheres = get_rotated_spheres_dict(in_sphere_file, _n_patches=args.n_patches)
    make_patch_fsaverage = make_patch_wrapper(spheres, _patch_size=args.patch_size,
                                              _resolution=args.resolution,
                                              method=args.interpolation)

    # We need to make sure 3 things line up before we zip them
    # 1) Special case: output is None, then put output with input
    if args.output_dirs == ['None']:
        args.output_dirs = args.input_dirs

    # 2) Input and output dirs
    li, lo, ls = len(args.input_dirs), len(args.output_dirs), len(args.scalar_files)
    if li > 1 and ls > 1 and not li == ls:
        raise ValueError(f'Input dir and scalar files could not be broadcast together with lengths ({li}, {ls})')
    elif lo > 1 and li > 1 and not lo == li:
        raise ValueError(f'Input dir and output dirs could not be broadcast together with lengths ({li}, {lo})')
    elif li > 1:
        # Prepare for iteration over input directories
        args.scalar_files = args.scalar_files * li
        if lo == 1:
            args.output_dirs = args.output_dirs * li
    elif ls > 1:
        # Prepare for iteration over scalar files
        args.input_dirs = args.input_dirs * ls
        if lo == 1:
            args.output_dirs = args.output_dirs * ls

    # 3) is_right parameter:
    if len(args.is_rights) < len(args.scalar_files):
        print(f'Warning: Not every scalar file was provided with hemisphere information.'
              f' Assuming they are all the same: {args.is_rights[0]}')
        args.is_rights = [args.is_rights[0]] * len(args.scalar_files)

    # For each scalar path that was supplied...
    for not_first, (scalar_file, input_dir, is_right, output_dir) in enumerate(zip(args.scalar_files,
                                                                                   args.input_dirs,
                                                                                   args.is_rights,
                                                                                   args.output_dirs)):
        makedirs(output_dir, exist_ok=True)
        # Create output names
        name_template = join(output_dir, scalar_file.replace('.gii', '').replace('.', '-'))
        if not args.stack_name:
            out_name_stack = name_template + '.npy'
            if isfile(out_name_stack) and not args.overwrite:
                print(f'FileExists: {out_name_stack} and --overwriting is False')
                continue
        else:
            out_name_stack = join(output_dir, args.stack_name)

        if not args.boxfg_name:
            out_name_boxfg = name_template + '.jpg'
        else:
            out_name_boxfg = join(output_dir, args.boxfg_name)

        out_name_prcfg = join(output_dir, args.prcfg_name)

        print(f'Computing patches for "{scalar_file}".')
        base, ext = splitext(scalar_file)
        if ext == '.gii':
            make_patch = make_patch_fsaverage
        else:
            # For non-fsaverage files (?h.thickness, ?h.curv) we read proprietary sphere data
            in_sphere_file = join(input_dir, base + '.sphere')
            spheres = get_rotated_spheres_dict(in_sphere_file, _n_patches=args.n_patches)
            make_patch = make_patch_wrapper(spheres, _patch_size=args.patch_size,
                                            _resolution=args.resolution,
                                            method=args.interpolation)
        in_scalar_file = join(input_dir, scalar_file)
        patch_stack, process_figs = make_stack(make_patch, _scalar_file=in_scalar_file,
                                               verbose=args.verbose and (not not_first),
                                               _is_right=is_right)

        # Visualization and storing (only do it for the first sample)
        if args.verbose and not not_first:
            print(f'Storing verbose figures to "{name_template}".')
        for j, proc_fig in enumerate(tqdm(process_figs)):
            proc_fig.suptitle(proc_fig._suptitle.get_text().format(name_template))
            proc_fig.savefig(out_name_prcfg.format(str(j).zfill(2)))

        # Create and store a box figure to aid quick interpretation of the stacked results
        box_fig = stack_fig(patch_stack, scalar_file, rotations=OrderedSet([rot for _, rot in spheres[0]]))
        box_fig.savefig(out_name_boxfg)
        plt.close(box_fig)
        np.save(out_name_stack, patch_stack)
        print(f'Successfully created: "{out_name_stack}" from "{in_sphere_file}"')
    sys.exit('sphere2dice.py complete successfully!')