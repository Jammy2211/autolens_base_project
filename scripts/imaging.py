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

    Dataset-specific properties are read from info.json, which lives alongside the dataset.
    See the README for the expected info.json format.
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

    Load dataset metadata from info.json. All fields use sensible defaults if absent,
    so real-data info.json files without every field will still work.

    See the README for a full description of the info.json format.
    """
    with open(dataset_path / "info.json") as json_file:
        info = json.load(json_file)

    pixel_scale = info.get("pixel_scale", 0.05)
    n_batch = info.get("n_batch", 50)
    mask_radius = info.get("mask_radius", 3.5)
    subhalo_grid_dimension_arcsec = info.get("subhalo_grid_dimensions_arcsec", 3.0)

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

    dataset = al.Imaging.from_fits(
        data_path=dataset_path / "data.fits",
        psf_path=dataset_path / "psf.fits",
        noise_map_path=dataset_path / "noise_map.fits",
        pixel_scales=pixel_scale,
    )

    positions = al.Grid2DIrregular(
        al.from_json(file_path=dataset_path / "positions.json")
    )

    mask = al.Mask2D.circular(
        shape_native=dataset.shape_native,
        pixel_scales=pixel_scale,
        centre=(0.0, 0.0),
        radius=mask_radius,
    )

    dataset = dataset.apply_mask(mask=mask)

    over_sample_size = al.util.over_sample.over_sample_size_via_radial_bins_from(
        grid=dataset.grid,
        sub_size_list=[4, 2, 1],
        radial_list=[0.1, 0.3],
    )

    dataset = dataset.apply_over_sampling(
        over_sample_size_lp=over_sample_size,
    )

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
    __SOURCE LP PIPELINE__
    """
    analysis = al.AnalysisImaging(
        dataset=dataset,
        positions_likelihood_list=[al.PositionsLH(threshold=0.4, positions=positions)],
    )

    lens_bulge = al.model_util.mge_model_from(
        mask_radius=mask_radius,
        total_gaussians=30,
        gaussian_per_basis=2,
        centre_prior_is_uniform=True,
    )

    source_bulge = al.model_util.mge_model_from(
        mask_radius=mask_radius, total_gaussians=20, centre_prior_is_uniform=False
    )

    source_lp_result = slam_pipeline.source_lp.run(
        settings_search=settings_search,
        analysis=analysis,
        lens_bulge=lens_bulge,
        lens_disk=None,
        mass=af.Model(al.mp.Isothermal),
        shear=af.Model(al.mp.ExternalShear),
        source_bulge=source_bulge,
        mass_centre=(0.0, 0.0),
        redshift_lens=redshift_lens,
        redshift_source=redshift_source,
    )

    """
    __SOURCE PIX PIPELINE 1__
    """
    hilbert_pixels = al.model_util.hilbert_pixels_from_pixel_scale(pixel_scale)

    image_mesh = al.image_mesh.Hilbert(pixels=hilbert_pixels, weight_power=3.5, weight_floor=0.01)

    galaxy_image_name_dict = al.galaxy_name_image_dict_via_result_from(
        result=source_lp_result
    )

    image_plane_mesh_grid = image_mesh.image_plane_mesh_grid_from(
        mask=dataset.mask, adapt_data=galaxy_image_name_dict["('galaxies', 'source')"]
    )

    edge_pixels_total = 30

    image_plane_mesh_grid = al.image_mesh.append_with_circle_edge_points(
        image_plane_mesh_grid=image_plane_mesh_grid,
        centre=mask.mask_centre,
        radius=mask_radius + mask.pixel_scale / 2.0,
        n_points=edge_pixels_total,
    )

    adapt_images = al.AdaptImages(
        galaxy_name_image_dict=galaxy_image_name_dict,
        galaxy_name_image_plane_mesh_grid_dict={
            "('galaxies', 'source')": image_plane_mesh_grid
        },
    )

    signal_to_noise_threshold = 3.0
    over_sample_size_pixelization = np.where(
        galaxy_image_name_dict["('galaxies', 'source')"] > signal_to_noise_threshold,
        4,
        2,
    )
    over_sample_size_pixelization = al.Array2D(
        values=over_sample_size_pixelization, mask=mask
    )

    dataset = dataset.apply_over_sampling(
        over_sample_size_lp=over_sample_size,
        over_sample_size_pixelization=over_sample_size_pixelization,
    )

    analysis = al.AnalysisImaging(
        dataset=dataset,
        adapt_images=adapt_images,
        positions_likelihood_list=[
            source_lp_result.positions_likelihood_from(
                factor=2.0, minimum_threshold=0.3, positions=positions
            )
        ],
    )

    source_pix_result_1 = slam_pipeline.source_pix.run_1(
        settings_search=settings_search,
        analysis=analysis,
        source_lp_result=source_lp_result,
        mesh_init=al.mesh.Delaunay(
            pixels=hilbert_pixels, zeroed_pixels=edge_pixels_total
        ),
        regularization_init=af.Model(al.reg.AdaptSplit),
        n_batch=n_batch,
    )

    """
    __SOURCE PIX PIPELINE 2__
    """
    galaxy_image_name_dict = al.galaxy_name_image_dict_via_result_from(
        result=source_pix_result_1
    )

    image_mesh = al.image_mesh.Hilbert(
        pixels=hilbert_pixels,
        weight_power=3.5,
        weight_floor=0.01,
    )

    signal_to_noise_threshold_image_mesh = 3.0
    adapt_data_snr_max = galaxy_image_name_dict["('galaxies', 'source')"]
    adapt_data_snr_max[adapt_data_snr_max > signal_to_noise_threshold_image_mesh] = (
        signal_to_noise_threshold_image_mesh
    )

    image_plane_mesh_grid = image_mesh.image_plane_mesh_grid_from(
        mask=dataset.mask, adapt_data=adapt_data_snr_max
    )

    image_plane_mesh_grid = al.image_mesh.append_with_circle_edge_points(
        image_plane_mesh_grid=image_plane_mesh_grid,
        centre=mask.mask_centre,
        radius=mask_radius + mask.pixel_scale / 2.0,
        n_points=edge_pixels_total,
    )

    adapt_images = al.AdaptImages(
        galaxy_name_image_dict=galaxy_image_name_dict,
        galaxy_name_image_plane_mesh_grid_dict={
            "('galaxies', 'source')": image_plane_mesh_grid
        },
    )

    over_sample_size_pixelization = np.where(
        galaxy_image_name_dict["('galaxies', 'source')"] > signal_to_noise_threshold,
        4,
        2,
    )
    over_sample_size_pixelization = al.Array2D(
        values=over_sample_size_pixelization, mask=mask
    )

    dataset = dataset.apply_over_sampling(
        over_sample_size_lp=over_sample_size,
        over_sample_size_pixelization=over_sample_size_pixelization,
    )

    analysis = al.AnalysisImaging(
        dataset=dataset,
        adapt_images=adapt_images,
    )

    source_pix_result_2 = slam_pipeline.source_pix.run_2(
        settings_search=settings_search,
        analysis=analysis,
        source_lp_result=source_lp_result,
        source_pix_result_1=source_pix_result_1,
        mesh=al.mesh.Delaunay(
            pixels=hilbert_pixels, zeroed_pixels=edge_pixels_total
        ),
        regularization=af.Model(al.reg.AdaptSplit),
        n_batch=n_batch,
    )

    """
    __LIGHT LP PIPELINE__
    """
    analysis = al.AnalysisImaging(
        dataset=dataset,
        adapt_images=adapt_images,
    )

    lens_bulge = al.model_util.mge_model_from(
        mask_radius=mask_radius,
        total_gaussians=30,
        gaussian_per_basis=2,
        centre_prior_is_uniform=True,
    )

    light_result = slam_pipeline.light_lp.run(
        settings_search=settings_search,
        analysis=analysis,
        source_result_for_lens=source_pix_result_1,
        source_result_for_source=source_pix_result_2,
        pixel_scales=pixel_scale,
        lens_bulge=lens_bulge,
        lens_disk=None,
        n_batch=n_batch,
    )

    """
    __MASS TOTAL PIPELINE__
    """
    analysis = al.AnalysisImaging(
        dataset=dataset,
        adapt_images=adapt_images,
        positions_likelihood_list=[
            source_pix_result_2.positions_likelihood_from(
                factor=3.0, minimum_threshold=0.2
            )
        ],
    )

    mass_result = slam_pipeline.mass_total.run(
        settings_search=settings_search,
        analysis=analysis,
        source_result_for_lens=source_pix_result_1,
        source_result_for_source=source_pix_result_2,
        light_result=light_result,
        mass=af.Model(al.mp.PowerLaw),
        n_batch=n_batch,
    )

    """
    __SUBHALO PIPELINE__
    """
    analysis = al.AnalysisImaging(
        dataset=dataset,
        adapt_images=adapt_images,
        positions_likelihood_list=[mass_result.positions_likelihood_from(
            factor=2.0, minimum_threshold=0.2,
        )],
    )

    subhalo_grid_search_result_1 = slam_pipeline.subhalo.detection.run_1_grid_search(
        settings_search=settings_search,
        analysis=analysis,
        mass_result=mass_result,
        subhalo_mass=af.Model(al.mp.NFWMCRLudlowSph),
        grid_dimension_arcsec=subhalo_grid_dimension_arcsec,
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

    parser = argparse.ArgumentParser(description="PyAutoLens SLAM Pipeline")

    parser.add_argument(
        "--sample",
        metavar="name",
        required=False,
        default=None,
        help="Name of the sample subdirectory inside dataset/ (e.g. slacs_sample).",
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
