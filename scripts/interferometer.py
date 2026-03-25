import numpy as np
import sys
import json
from pathlib import Path


def fit(dataset_name, sample_name):
    """
    __Paths__

    The project root is inferred from the location of this script, so no paths need to be
    hard-coded. This means the script runs identically whether invoked locally or on the HPC —
    the batch submission script simply calls this file with its absolute path.

    Dataset-specific properties (pixel_scale, n_batch) are read from info.json, which lives
    alongside the dataset. See the README for the expected info.json format.
    """
    project_root = Path(__file__).parent.parent

    config_path = project_root / "config"
    dataset_path = project_root / "dataset"
    output_path = project_root / "output"

    if sample_name is not None:
        dataset_path = dataset_path / sample_name
    if dataset_name is not None:
        dataset_path = dataset_path / dataset_name

    """
    __Info__

    Load dataset metadata from info.json. This file must contain at minimum:

        pixel_scale       -- arcseconds per pixel (float)
        n_batch           -- pixelization batch size, lower for higher-resolution data (int)

    Optional fields:

        real_space_shape  -- [height, width] of the real-space mask (default [256, 256])
        mask_radius       -- circular mask radius in arcseconds (default 3.5)

    See the README for a full description of the info.json format.
    """
    with open(dataset_path / "info.json") as json_file:
        info = json.load(json_file)

    pixel_scale = info["pixel_scale"]
    n_batch = info["n_batch"]
    real_space_shape = info.get("real_space_shape", [256, 256])
    mask_radius = info.get("mask_radius", 3.5)

    """
    __Configs__
    """
    from autoconf import conf

    conf.instance.push(new_path=config_path, output_path=output_path)

    """
    __AutoLens + Data__
    """
    import autofit as af
    import autolens as al

    sys.path.insert(0, str(project_root))
    import slam_pipeline

    real_space_mask = al.Mask2D.circular(
        shape_native=tuple(real_space_shape),
        pixel_scales=pixel_scale,
        radius=mask_radius,
    )

    dataset = al.Interferometer.from_fits(
        data_path=dataset_path / "data.fits",
        noise_map_path=dataset_path / "noise_map.fits",
        uv_wavelengths_path=dataset_path / "uv_wavelengths.fits",
        real_space_mask=real_space_mask,
        transformer_class=al.TransformerNUFFT,
    )

    positions = al.Grid2DIrregular(
        al.from_json(file_path=dataset_path / "positions.json")
    )

    try:
        nufft_precision_operator = np.load(
            file=dataset_path / "nufft_precision_operator.npy",
        )
    except FileNotFoundError:
        nufft_precision_operator = None

    dataset = dataset.apply_sparse_operator(
        nufft_precision_operator=nufft_precision_operator, use_jax=True, show_progress=True
    )

    """
    __Positions Likelihood__
    """
    positions_likelihood = al.PositionsLH(positions=positions, threshold=0.3)

    """
    __Settings__

    Disable the positive-only linear algebra solver so the source reconstruction can
    have negative pixel values, which is required for interferometer data.
    """
    settings = al.Settings(use_positive_only_solver=False)

    """
    __Settings AutoFit__
    """
    settings_search = af.SettingsSearch(
        path_prefix=sample_name,
        unique_tag=dataset_name,
        session=None,
        info=info,
    )

    """
    __Redshifts__
    """
    redshift_lens = info.get("redshift_lens", 0.5)
    redshift_source = info.get("redshift_source", 1.0)

    """
    __SOURCE PIX PIPELINE 1__

    Uses an Overlay image mesh to build an image-plane mesh grid that covers the
    entire real-space mask. An Overlay grid is used (rather than Hilbert) because
    no adapt image is available at this stage — there is no SOURCE LP stage for
    interferometer data.
    """
    image_mesh = al.image_mesh.Overlay(shape=(26, 26))

    image_plane_mesh_grid = image_mesh.image_plane_mesh_grid_from(
        mask=dataset.mask,
    )

    edge_pixels_total = 30

    image_plane_mesh_grid = al.image_mesh.append_with_circle_edge_points(
        image_plane_mesh_grid=image_plane_mesh_grid,
        centre=real_space_mask.mask_centre,
        radius=mask_radius + real_space_mask.pixel_scale / 2.0,
        n_points=edge_pixels_total,
    )

    adapt_images = al.AdaptImages(
        galaxy_name_image_plane_mesh_grid_dict={
            "('galaxies', 'source')": image_plane_mesh_grid
        },
    )

    analysis = al.AnalysisInterferometer(
        dataset=dataset,
        adapt_images=adapt_images,
        positions_likelihood_list=[positions_likelihood],
        settings=settings,
    )

    source_pix_result_1 = slam_pipeline.source_pix.run_1__bypass_lp(
        settings_search=settings_search,
        analysis=analysis,
        lens_bulge=None,
        lens_disk=None,
        mass=af.Model(al.mp.Isothermal),
        shear=af.Model(al.mp.ExternalShear),
        mesh_init=al.mesh.Delaunay(
            pixels=hilbert_pixels, zeroed_pixels=edge_pixels_total
        ),
        regularization_init=af.Model(al.reg.ConstantSplit),
        redshift_lens=redshift_lens,
        redshift_source=redshift_source,
        n_batch=n_batch,
    )

    """
    __SOURCE PIX PIPELINE 2__

    The lens light pipeline is skipped for interferometer data because the lens light
    is not present in visibility data.
    """
    galaxy_image_name_dict = al.galaxy_name_image_dict_via_result_from(
        result=source_pix_result_1, use_model_images=True
    )

    image_mesh = al.image_mesh.Hilbert(pixels=hilbert_pixels, weight_power=3.5, weight_floor=0.01)

    image_plane_mesh_grid = image_mesh.image_plane_mesh_grid_from(
        mask=dataset.mask, adapt_data=galaxy_image_name_dict["('galaxies', 'source')"]
    )

    image_plane_mesh_grid = al.image_mesh.append_with_circle_edge_points(
        image_plane_mesh_grid=image_plane_mesh_grid,
        centre=real_space_mask.mask_centre,
        radius=mask_radius + real_space_mask.pixel_scale / 2.0,
        n_points=edge_pixels_total,
    )

    adapt_images = al.AdaptImages(
        galaxy_name_image_dict=galaxy_image_name_dict,
        galaxy_name_image_plane_mesh_grid_dict={
            "('galaxies', 'source')": image_plane_mesh_grid
        },
    )

    analysis = al.AnalysisInterferometer(
        dataset=dataset,
        adapt_images=adapt_images,
        settings=settings,
    )

    source_pix_result_2 = slam_pipeline.source_pix.run_2(
        settings_search=settings_search,
        analysis=analysis,
        source_lp_result=source_pix_result_1,
        source_pix_result_1=source_pix_result_1,
        mesh=al.mesh.Delaunay(
            pixels=hilbert_pixels, zeroed_pixels=edge_pixels_total
        ),
        regularization=af.Model(al.reg.AdaptSplit),
        n_batch=n_batch,
    )

    """
    __MASS TOTAL PIPELINE__
    """
    analysis = al.AnalysisInterferometer(
        dataset=dataset,
        adapt_images=adapt_images,
        positions_likelihood_list=[
            source_pix_result_1.positions_likelihood_from(
                factor=3.0, minimum_threshold=0.2
            )
        ],
        settings=settings,
    )

    mass_result = slam_pipeline.mass_total.run(
        settings_search=settings_search,
        analysis=analysis,
        source_result_for_lens=source_pix_result_1,
        source_result_for_source=source_pix_result_2,
        light_result=None,
        mass=af.Model(al.mp.PowerLaw),
        multipole_1=af.Model(al.mp.PowerLawMultipole),
        multipole_3=af.Model(al.mp.PowerLawMultipole),
        multipole_4=af.Model(al.mp.PowerLawMultipole),
        n_batch=n_batch,
    )

    """
    __SUBHALO PIPELINE__
    """
    analysis = al.AnalysisInterferometer(
        dataset=dataset,
        adapt_images=adapt_images,
        positions_likelihood_list=[
            mass_result.positions_likelihood_from(
                factor=3.0, minimum_threshold=0.2
            )
        ],
        settings=settings,
    )

    subhalo_grid_search_result_1 = slam_pipeline.subhalo.detection.run_1_grid_search(
        settings_search=settings_search,
        analysis=analysis,
        mass_result=mass_result,
        subhalo_mass=af.Model(al.mp.NFWMCRLudlowSph),
        grid_dimension_arcsec=3.0,
        number_of_steps=2,
    )

    slam_pipeline.subhalo.detection.run_2_subhalo(
        settings_search=settings_search,
        analysis=analysis,
        mass_result=mass_result,
        subhalo_grid_search_result_1=subhalo_grid_search_result_1,
        subhalo_mass=af.Model(al.mp.NFWMCRLudlowSph),
    )


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="PyAutoLens SLAM Interferometer Pipeline")

    parser.add_argument(
        "--sample",
        metavar="name",
        required=False,
        default=None,
        help="Name of the sample subdirectory inside dataset/ (e.g. alma_sample).",
    )

    parser.add_argument(
        "--dataset",
        metavar="name",
        required=False,
        default=None,
        help="Name of the dataset subdirectory inside dataset/[sample]/ (e.g. lens0001).",
    )

    args = parser.parse_args()

    fit(dataset_name=args.dataset, sample_name=args.sample)

"""
Finish.
"""
