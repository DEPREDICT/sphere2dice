# What does this code do?

It converts Cortical Brain Morphometry data from FreeSurfer (`lh.thickness`) or Gifty (`lh-thickness.gii`) type into a stack of 2D images and saves is as numpy array. This stack of images is easier to use as input to conventional CNNs, while avoiding projection issues other methods suffer from.
![repository-open-graph-template](https://user-images.githubusercontent.com/38399483/158713220-cde45111-8d71-4e09-8c55-6516e6749ba8.png)
# How to use?

Run `python sphere2dice.py -h` for help on the inputs.

To get started processing one single file of your own you only have to provide the `--scalar_file` parameter.
`python sphere2dice.py --scalar_files sub-01_Thickness.gii `

This of course also works in series:

`python sphere2dice.py --scalar_files sub-01_Thickness.gii sub-02_Thickness.gii` 

Alternatively, it is also possible to iterate over directories:

`python sphere2dice.py --scalar_files lh.thickness --input_dirs sub-001/Surf sub-002/Surf-sub-003/Surf output_dirs sub-001 sub-002 sub-003`

If we want to store our results in the same folder as we read them from, we can set `output_dirs` to `None`:

`python sphere2dice.py --scalar_files lh.thickness --input_dirs sub-001/Surf sub-002/Surf-sub-003/Surf output_dirs None`

# How sphere2dice works
It uses fsaverage sphere coordinates (provided in separate file as `lh.sphere`) for fsaverage registered Gifty (.gii) files. Alternatively, for FreeSurfer Morphometry data, it is assumed that the scalar file read is accompanied by a corresponding `?h.sphere` file.

The index of these sphere coordinates is linked to the index of scalar values (like thickness or curvature).
These coordinates are projected to a 2D plane, and we cut a square out of the center of size `size patch_size`.
Then we interpolate these points to a given `resolution`. 
We create multiple patches from all different sides but rotating the coordinates and repeating the process.
If the code is run with the `--verbose` flag, extra figures are generated displaying this process.

# Requirements
Compatible with Python 3. 
Modules used are basic so requirements are minimal.
I have therefore not containerized the script.
I did attach the `requirements.txt` for repeatability. 

# Terminology
* *Patch*: A 2D projection. You have two choices for the number of patches as of now: 6 or 18. This makes sense as 6 patches contain two opposite projections from three orthogonal axes, like folding a simple box or dice. 
  An issue with this lies in the edges: As when folding a 6-sided die: there is no overlap and the anatomy at the edges lacks anatomical context. The solution is simple: Take an extra patch right in between every two patches, this comes down to 4 around the equator, 4 on the north pole and 4 on the south pole = 18.
  The code contains debugging modes that show this process (but might be a bit finicky to set up, so I’ll supply these as soon as I’m back at my machine.
* *Dice*: I use the word dice every so often (like in the name of the code) for folding a simple box, which is the most simple way to perceive this method (but it is not limited to it)
* *Stack*: A stack is a 3D volume of 2D patches. 

# Example code

## Run the demo
Create additional figures using the `verbose` flag displaying the rotation process.
See `/output/angle-0*_overview.jpg` for the results.

    (ENIGMA) D:\...\sphere2dice-standalone>python sphere2dice.py --verbose
    Computing patches for "example-thickness-file.gii".
    100%|████████████████████████████████████████████████████████████████████████████| 6/6 [00:13<00:00,  2.28s/it]
    Storing verbose figures to "out_name_prcfg".
    100%|████████████████████████████████████████████████████████████████████████████| 6/6 [01:36<00:00, 16.02s/it]
    Successfully created: "output\example-thickness-file.npy" from "input\lh.sphere"
    sphere2dice.py complete succesfully

## Real world example

There are two usefull ways to use the script:

*1) Iterate over files in the same folder*. Note the use of the `is_right` flag.

python shere2dice.py --scalar_files lh.thickness --input_dirs sub-001 sub-002 --output_dirs None

    (ENIGMA) D:\...\sphere2dice-standalone>python sphere2dice.py --scalar_files 26567_space-fsaverage_hemi-R.thickness.gii 26677_space-fsaverage_hemi-R.thickness.gii --input_dir D:\...\AnatSurfRH\Thickness --is_right True True
    Computing patches for "26567_space-fsaverage_hemi-R.thickness.gii".
    100%|████████████████████████████████████████████████████████████████████████████| 6/6 [00:11<00:00,  1.88s/it]
    0it [00:00, ?it/s]
    Successfully created: "output\26567_space-fsaverage_hemi-R.thickness.npy" from "input\lh.sphere"
    Computing patches for "26677_space-fsaverage_hemi-R.thickness.gii".
    100%|████████████████████████████████████████████████████████████████████████████| 6/6 [00:11<00:00,  1.93s/it]
    0it [00:00, ?it/s]
    Successfully created: "output\26677_space-fsaverage_hemi-R.thickness.npy" from "input\lh.sphere"
    sphere2dice.py complete succesfully

*2) Iterate over folders*: Same file (`lh.thickness`) in different input folders (`sub-001`,  `sub-002`), with `--output_dirs` set to `None`  output is stored in the respective input folders

```
python shere2dice.py --scalar_files lh.thickness --input_dirs sub-001 sub-002 --output_dirs None
```

# Origin of accompanying files
Origin of the lh.sphere[^1] file: CorticalParcellation_Yeo2011 - Free Surfer Wiki [^2]
Then `Yeo_JNeurophysiol11_FreeSurfer\fsaverage\surf\lh.sphere`

[^1]: ftp://surfer.nmr.mgh.harvard.edu/pub/data/Yeo_JNeurophysiol11_FreeSurfer.zip
[^2]: https://surfer.nmr.mgh.harvard.edu/fswiki/CorticalParcellation_Yeo2011
