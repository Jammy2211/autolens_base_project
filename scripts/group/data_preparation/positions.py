"""
GUI Preprocessing: Positions
============================

This tool allows one to input the positions of strong lenses via a GUI, which can be used to resample inaccurate
mass models during lensing modeling.

This GUI is adapted from the following code: https://gist.github.com/brikeats/4f63f867fd8ea0f196c78e9b835150ab
"""
# %matplotlib inline
# from pyprojroot import here
# workspace_path = str(here())
# %cd $workspace_path
# print(f"Working Directory has been set to `{workspace_path}`")

import os
from os import path
import autolens as al
import autolens.plot as aplt
from matplotlib import pyplot as plt
from scipy.ndimage import center_of_mass
import numpy as np

#from autolens_workspace.scripts.plot.visuals_2d.PositionsScatter import positions

"""
__Dataset__

Setup the path the datasets we'll use to illustrate preprocessing, which is the 
folder `dataset/imaging/simple__no_lens_light`.
"""
dataset_name = "102021990_NEG650312660474055399"
dataset_path = path.join("..","..","..", "dataset", "sample_group", dataset_name)#, "F444W")

"""
The pixel scale of the imaging dataset.
"""
pixel_scales = 0.1

"""
Load the image which we will use to mark the positions.
"""
data = al.Array2D.from_fits(
    file_path=path.join(dataset_path, "data.fits"), pixel_scales=pixel_scales#, hdu=1,
)

# vis_index = 0
#
# data = al.Imaging.from_fits(
#     data_path=path.join(dataset_path, "data.fits"),
#     data_hdu=vis_index + 1,
#     noise_map_path=path.join(dataset_path, "data.fits"),
#     noise_map_hdu=vis_index + 3,
#     psf_path=path.join(dataset_path, "data.fits"),
#     psf_hdu=vis_index + 2,
#     pixel_scales=0.1,
#     check_noise_map=False,
# )


mask_radius = 6.5

mask = al.Mask2D.circular(
    shape_native=data.shape_native,
    pixel_scales=data.pixel_scales,
    radius=mask_radius,
)

data = data.apply_mask(mask=mask)

# lens_light = al.Array2D.from_fits(
#      file_path=path.join(dataset_path, "result", "mge_lens_light.fits"), pixel_scales=pixel_scales
# )

# lens_light = al.Array2D.from_fits(
#     file_path=path.join(dataset_path, "model_lens_light.fits"), pixel_scales=pixel_scales
# )

# data = data - lens_light

# Apply log10 transformation while ensuring no zero or negative values
log_data = np.log10(np.maximum(data.native, 1e-6))

# Convert back to an Array2D object
log_data = al.Array2D.no_mask(values=log_data, pixel_scales=pixel_scales)

"""
__Search Box__

When you click on a pixel to mark a position, the search box looks around this click and finds the pixel with
the highest flux to mark the position.

The `search_box_size` is the number of pixels around your click this search takes place.
"""
search_box_size = 2

"""
__Clicker__

Set up the `Clicker` object from the `clicker.py` module, which monitors your mouse clicks in order to determine
the positions.
"""
clicker = al.Clicker(
    image=log_data, pixel_scales=pixel_scales, search_box_size=search_box_size
)


"""
For lenses with bright lens light emission, it can be difficult to get the source light to show. The normalization
below uses a log-scale with a capped maximum, which better contrasts the lens and source emission.
"""
cmap = aplt.Cmap(norm="linear", vmin=1.0e-4, vmax=np.max(data))

norm = cmap.norm_from(array=None)

"""
Set up the clicker canvas and load the GUI which you can now click on to mark the positionss.
"""
n_y, n_x = data.shape_native
hw = n_x / 2 * pixel_scales
ext = [-hw, hw, -hw, hw]
fig = plt.figure(figsize=(14, 14))
plt.imshow(log_data.native, cmap="jet", extent=ext)
plt.colorbar()
cid = fig.canvas.mpl_connect("button_press_event", clicker.onclick)
plt.show()
fig.canvas.mpl_disconnect(cid)
plt.close(fig)


"""
Use the results of the Clicker GUI to create the list of the positions.
"""
positions_pixel_centres = al.Grid2DIrregular(values=clicker.click_list)

def physical_to_pixel(physical_coordinates, pixel_size, image_size, image_extent):
    """
    Convert physical coordinates to pixel coordinates.

    Parameters:
    - physical_coordinates: (y, x) tuple in physical units (e.g., arcseconds).
    - pixel_size: Size of each pixel in physical units (e.g., arcseconds/pixel).
    - image_size: Tuple (height, width) representing the image size in pixels.
    - image_extent: Tuple (min_y, max_y, min_x, max_x) representing the physical extent of the image.

    Returns:
    - Tuple (y_pixel, x_pixel) of pixel coordinates.
    """
    # Get the physical coordinate bounds of the image
    min_y, max_y, min_x, max_x = image_extent

    # Convert physical coordinates to pixel coordinates
    y_pixel = int((physical_coordinates[0] - min_y) / pixel_size)
    x_pixel = int((physical_coordinates[1] - min_x) / pixel_size)

    # Ensure the pixel coordinates are within the image bounds
    y_pixel = np.clip(y_pixel, 0, image_size[0] - 1)
    x_pixel = np.clip(x_pixel, 0, image_size[1] - 1)

    return y_pixel, x_pixel


def find_subpixel_centroid(image, center_physical, pixel_size, image_extent, window_size=3):
    """
    Find the sub-pixel centroid of the brightest region around a given physical location.

    Parameters:
    - image: 2D numpy array representing the image data (magnification map).
    - center_physical: Tuple (y, x) physical coordinates of the brightest region.
    - pixel_size: Size of each pixel in physical units (e.g., arcseconds/pixel).
    - image_extent: Tuple (min_y, max_y, min_x, max_x) representing the physical extent of the image.
    - window_size: Size of the square window around the initial brightest region to consider for centroiding.

    Returns:
    - Sub-pixel accurate (y, x) centroid in physical coordinates.
    """
    # Convert physical coordinates to pixel coordinates
    center_pixel = physical_to_pixel(center_physical, pixel_size, image.shape_native, image_extent)

    # Define the window bounds around the center pixel
    half_window = window_size // 2
    y_min = max(center_pixel[0] - half_window, 0)
    y_max = min(center_pixel[0] + half_window + 1, image.shape_native[0])
    x_min = max(center_pixel[1] - half_window, 0)
    x_max = min(center_pixel[1] + half_window + 1, image.shape_native[1])

    # Extract the image region around the center pixel
    sub_image = image.native[::-1, :][y_min:y_max, x_min:x_max] # y-axis has to be flipped

    # Find the centroid of the region using the center of mass
    # This gives sub-pixel precision since it computes the weighted center of the image region
    y_centroid, x_centroid = center_of_mass(sub_image)

    # Convert the centroid back to physical coordinates
    y_physical = (y_centroid + 0.5 + y_min) * pixel_size + image_extent[0] # The sub_image origin is placed at the centre of the first pixel, not the left-lower corner, hence plus 0.5
    x_physical = (x_centroid + 0.5 + x_min) * pixel_size + image_extent[2]

    return y_physical, x_physical

subpixel_centroid = []
for position_pixel_centre in positions_pixel_centres:
    subpixel_centroid.append(find_subpixel_centroid(data, position_pixel_centre, pixel_scales, ext, window_size=3))

positions = al.Grid2DIrregular(values=subpixel_centroid)

"""
__Output__

Now lets plot the image and positions,, so we can check that the positions overlap the brightest pixels in the
lensed source.
"""
visuals = aplt.Visuals2D(mass_profile_centres=(positions, positions_pixel_centres))

array_2d_plotter = aplt.Array2DPlotter(
    array=log_data, visuals_2d=visuals, mat_plot_2d=aplt.MatPlot2D()
)
array_2d_plotter.figure_2d()

"""
Output this image of the positions to a .png file in the dataset folder for future reference.
"""
array_2d_plotter = aplt.Array2DPlotter(
    array=log_data,
    visuals_2d=visuals,
    mat_plot_2d=aplt.MatPlot2D(
        output=aplt.Output(path=dataset_path, filename="positions", format="png")
    ),
)
array_2d_plotter.figure_2d()

"""
Output the positions to a .json file in the dataset folder, so we can load them in modeling scripts.
"""
al.output_to_json(
    obj=positions,
    file_path=path.join(dataset_path, "positions.json"),
)
