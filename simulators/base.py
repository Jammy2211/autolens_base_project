import numpy as np
import json
import sys
from pathlib import Path

import autolens as al
import autolens.plot as aplt

"""
__Paths__

The project root is inferred from the location of this script, so no paths need to be
hard-coded. This means the script runs identically whether invoked locally or on the HPC.
"""
project_root = Path(__file__).parent.parent

config_path = project_root / "config"
dataset_path = project_root / "dataset"

dataset_name = sys.argv[1] if len(sys.argv) > 1 else None

if dataset_name is not None:
    dataset_path = dataset_path / dataset_name

"""
__Dataset Properties__

These define the properties of the simulated dataset. After simulation they are written to
info.json alongside the data, so that analysis scripts can read them back without any
hard-coding or duplication.

    pixel_scale  -- arcseconds per pixel
    shape_native -- image dimensions in pixels
    n_batch      -- pixelization batch size; lower values reduce VRAM for higher-resolution data
"""
pixel_scale = 0.05
shape_native = (400, 400)
n_batch = 40

"""
__Simulate__
"""
psf = al.Convolver.from_gaussian(
    shape_native=(21, 21), sigma=pixel_scale, pixel_scales=pixel_scale, normalize=True
)

grid = al.Grid2D.uniform(
    shape_native=shape_native,
    pixel_scales=pixel_scale,
)

simulator = al.SimulatorImaging(
    exposure_time=300.0,
    psf=psf,
    use_real_space_convolution=True,
    background_sky_level=0.1,
    add_poisson_noise_to_data=True,
)

source_galaxy = al.Galaxy(
    redshift=1.0,
    bulge=al.lp.SersicCore(
        centre=(0.01, 0.01),
        ell_comps=(0.096225, -0.055555),
        intensity=1.5,
        effective_radius=0.3,
        sersic_index=2.5,
    ),
)

lens_galaxy = al.Galaxy(
    redshift=0.5,
    bulge=al.lp.Sersic(
        centre=(0.0, 0.0),
        ell_comps=(0.0, -0.05),
        intensity=0.35,
        effective_radius=0.7,
        sersic_index=1.2,
    ),
    disk=al.lp.Sersic(
        centre=(0.0, 0.0),
        ell_comps=(0.02, -0.05),
        intensity=0.35,
        effective_radius=4.5,
        sersic_index=1.2,
    ),
    mass=al.mp.PowerLaw(
        centre=(0.0, 0.0),
        einstein_radius=2.0,
        ell_comps=(0.02, -0.02),
        slope=2.0,
    ),
    shear=al.mp.ExternalShear(gamma_1=0.0, gamma_2=0.05),
)

tracer = al.Tracer(galaxies=[lens_galaxy, source_galaxy])

"""
__Over Sampling__
"""
over_sample_size_lens = al.util.over_sample.over_sample_size_via_radial_bins_from(
    grid=grid,
    sub_size_list=[32, 8, 4],
    radial_list=[0.3, 0.6],
    centre_list=[(0.0, 0.0)],
)

traced_grid = tracer.traced_grid_2d_list_from(grid=grid)[-1]

over_sample_size_source = al.util.over_sample.over_sample_size_via_radial_bins_from(
    grid=traced_grid,
    sub_size_list=[32, 16, 4],
    radial_list=[0.3, 0.6],
    centre_list=[source_galaxy.bulge.centre],
)

over_sample_size = np.where(
    over_sample_size_source > over_sample_size_lens,
    over_sample_size_source,
    over_sample_size_lens,
)
over_sample_size = al.Array2D(values=over_sample_size, mask=grid.mask)

grid = grid.apply_over_sampling(over_sample_size=over_sample_size)

"""
__Output__
"""
dataset_path.mkdir(parents=True, exist_ok=True)

mat_plot_2d = aplt.MatPlot2D(
    output=aplt.Output(path=dataset_path, format=["png", "fits"])
)

tracer_plotter = aplt.Tracer(tracer=tracer, grid=grid, mat_plot_2d=mat_plot_2d)
tracer_plotter.figures_2d(image=True)

dataset = simulator.via_tracer_from(tracer=tracer, grid=grid)

mat_plot_2d = aplt.MatPlot2D(output=aplt.Output(path=dataset_path, format="png"))
imaging_plotter = aplt.Imaging(dataset=dataset, mat_plot_2d=mat_plot_2d)
imaging_plotter.subplot_dataset()

dataset.output_to_fits(
    data_path=dataset_path / "data.fits",
    psf_path=dataset_path / "psf.fits",
    noise_map_path=dataset_path / "noise_map.fits",
    overwrite=True,
)

al.output_to_json(obj=tracer, file_path=dataset_path / "tracer.json")

"""
__No-Lens Image__
"""
lens_image = lens_galaxy.padded_image_2d_from(
    grid=grid,
    psf_shape_2d=psf.shape_native,
)

lens_image = psf.convolved_image_from(image=lens_image, blurring_image=None)
lens_image = lens_image.trimmed_after_convolution_from(kernel_shape=psf.shape_native)

data_no_lens = dataset.data - lens_image

plotter = aplt.Array2DPlotter(array=data_no_lens, mat_plot_2d=mat_plot_2d)
plotter.set_filename("data_no_lens")
plotter.figure_2d()

data_no_lens.output_to_fits(
    file_path=dataset_path / "data_no_lens.fits", overwrite=True
)

"""
__Positions__
"""
solver = al.PointSolver.for_grid(
    grid=grid, pixel_scale_precision=0.001, magnification_threshold=0.1
)

positions = solver.solve(
    tracer=tracer, source_plane_coordinate=source_galaxy.bulge.centre
)

al.output_to_json(
    file_path=dataset_path / "positions.json",
    obj=positions,
)

"""
__info.json__

Write dataset properties to info.json. These are read by analysis scripts so that
pixel_scale and n_batch do not need to be hard-coded or passed as arguments — they
live with the data that defines them.
"""
info = {
    "pixel_scale": pixel_scale,
    "shape_native": list(shape_native),
    "n_batch": n_batch,
    "redshift_lens": lens_galaxy.redshift,
    "redshift_source": source_galaxy.redshift,
}

with open(dataset_path / "info.json", "w") as f:
    json.dump(info, f, indent=4)

print(f"Simulated dataset written to {dataset_path}")
