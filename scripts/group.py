"""
SLaM (Source, Light and Mass): Group
=====================================

This script models group-scale strong lenses that include multiple extra galaxies
whose light and mass contribute to the lens model. It uses the same four-stage SLaM
structure as `imaging.py` but extends it to handle two populations of extra galaxies:

- **Inner extra galaxies**: modeled with individual Isothermal mass profiles whose
  Einstein radii are bounded by their fitted luminosities.
- **Outer extra galaxies**: modeled with a shared luminosity scaling relation, where
  the Einstein radius is parameterized as ``scaling_factor * luminosity^scaling_exponent``.

The lens centre is determined automatically from the brightest pixel within a small
central region of the image, and an optional noise-scaling mask suppresses contaminating
flux from the extra galaxies during the initial light fit.

__This Script__

 - Lens light: MGE basis (30 Gaussians, 2 bases).
 - Lens mass: PowerLaw + ExternalShear + multipoles.
 - Source: Pixelization (Delaunay mesh, AdaptSplit regularization).
 - Extra galaxies: MGE light + Isothermal mass (inner) or scaling-relation mass (outer).

Pipelines: source_lp (run_0 + run_1) → source_pix → light_lp → mass_total.
"""
import numpy as np
import sys
import json
from pathlib import Path


def fit(dataset_name, sample_name):
    """
    __Paths__

    The project root is inferred from the location of this script, so no paths need to
    be hard-coded. This means the script runs identically whether invoked locally or on
    the HPC — the batch submission script simply calls this file with its absolute path.

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

    pixel_scale = info.get("pixel_scale", 0.1)
    mask_radius = info.get("mask_radius", 6.0)
    mask_centre = info.get("mask_centre", (0.0, 0.0))
    redshift_lens = info.get("redshift_lens", 0.5)
    redshift_source = info.get("redshift_source", 1.0)

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

    dataset_centre = dataset.data.brightest_sub_pixel_coordinate_in_region_from(
        region=(-0.3, 0.3, -0.3, 0.3), box_size=2
    )

    try:
        mask_extra_galaxies = al.Mask2D.from_fits(
            file_path=dataset_path / "mask_extra_galaxies.fits",
            pixel_scales=pixel_scale,
            invert=True,
        )
        dataset = dataset.apply_noise_scaling(mask=mask_extra_galaxies)
    except FileNotFoundError:
        pass

    """
    __Galaxies Centres and Positions__

    Inner extra galaxies are those close to the main lens; outer extra galaxies are
    further away and modeled via a shared luminosity scaling relation. Both populations
    are combined into a single list for over-sampling and masking purposes.
    """
    try:
        inner_extra_galaxies_centres = al.Grid2DIrregular(
            al.from_json(file_path=dataset_path / "inner_extra_galaxies_centres.json")
        )
    except FileNotFoundError:
        inner_extra_galaxies_centres = al.Grid2DIrregular([])

    try:
        outer_extra_galaxies_centres = al.Grid2DIrregular(
            al.from_json(file_path=dataset_path / "outer_extra_galaxies_centres.json")
        )
    except FileNotFoundError:
        outer_extra_galaxies_centres = al.Grid2DIrregular([])

    extra_galaxies_centres = al.Grid2DIrregular(
        inner_extra_galaxies_centres.in_list + outer_extra_galaxies_centres.in_list
    )

    positions = al.Grid2DIrregular(
        al.from_json(file_path=dataset_path / "positions.json")
    )

    lens_centre = dataset_centre

    """
    __Mask__
    """
    centre_list = [lens_centre] + list(extra_galaxies_centres)

    mask = al.Mask2D.circular(
        shape_native=dataset.shape_native,
        pixel_scales=dataset.pixel_scales,
        radius=mask_radius,
        centre=mask_centre,
    )

    dataset = dataset.apply_mask(mask=mask)

    over_sample_size = al.util.over_sample.over_sample_size_via_radial_bins_from(
        grid=dataset.grid,
        sub_size_list=[4, 2, 1],
        radial_list=[0.1, 0.3],
        centre_list=centre_list,
    )

    dataset = dataset.apply_over_sampling(over_sample_size_lp=over_sample_size)

    """
    __Settings AutoFit__
    """
    settings_search = af.SettingsSearch(
        path_prefix=sample_name,
        unique_tag=dataset_name,
        info=info,
        session=None,
    )


    """
    __SOURCE LP PIPELINE__

    The SOURCE LP PIPELINE runs in two stages.

    The first search (run_0) fits only the lens and extra galaxy light profiles with no
    mass model. This initializes the light model and measures the extra galaxy luminosities
    before the mass model is introduced.

    The second search (run_1) introduces the lens mass and source light, holding the extra
    galaxy light profiles fixed to their run_0 values. Inner extra galaxy mass profiles are
    initialized with Einstein radii bounded by their run_0 luminosities; outer extra galaxy
    mass profiles use a shared luminosity scaling relation.
    """
    analysis = al.AnalysisImaging(dataset=dataset)

    lens_bulge = al.model_util.mge_model_from(
        mask_radius=mask_radius,
        total_gaussians=30,
        gaussian_per_basis=2,
        centre_prior_is_uniform=False,
        centre=lens_centre,
        centre_sigma=0.1,
    )

    extra_galaxies_list = []

    for extra_galaxy_centre in extra_galaxies_centres:
        extra_galaxy_bulge = al.model_util.mge_model_from(
            mask_radius=mask_radius,
            total_gaussians=10,
            centre_prior_is_uniform=False,
            centre=(extra_galaxy_centre[0], extra_galaxy_centre[1]),
            centre_sigma=0.1,
        )
        extra_galaxy = af.Model(
            al.Galaxy, redshift=redshift_lens, bulge=extra_galaxy_bulge,
        )
        extra_galaxies_list.append(extra_galaxy)

    extra_galaxies_free_centres = af.Collection(extra_galaxies_list)

    source_lp_result_0 = slam_pipeline.source_lp.run_0__group(
        settings_search=settings_search,
        analysis=analysis,
        lens_bulge=lens_bulge,
        redshift_lens=redshift_lens,
        extra_galaxies=extra_galaxies_free_centres,
    )

    positions_likelihood = al.PositionsLH(positions=positions, threshold=0.3)

    analysis = al.AnalysisImaging(
        dataset=dataset, positions_likelihood_list=[positions_likelihood]
    )

    lens_mass = af.Model(al.mp.Isothermal)
    lens_mass.einstein_radius = af.UniformPrior(lower_limit=0.0, upper_limit=5.0)

    source_bulge = al.model_util.mge_model_from(
        mask_radius=1.0,
        total_gaussians=30,
        centre_prior_is_uniform=False,
        centre=lens_centre,
        centre_sigma=0.6,
    )

    extra_galaxies_fixed_list = []

    scaling_factor = af.UniformPrior(lower_limit=0.0, upper_limit=0.5)
    scaling_relation = af.UniformPrior(lower_limit=0.0, upper_limit=2.0)

    for i in range(len(inner_extra_galaxies_centres)):
        extra_galaxy_bulge = source_lp_result_0.instance.extra_galaxies[i].bulge

        mass = af.Model(al.mp.Isothermal)
        mass.centre = extra_galaxy_bulge.centre
        mass.ell_comps = extra_galaxy_bulge.ell_comps

        luminosity_per_gaussian_list = [
            2 * np.pi * g.sigma ** 2 / g.axis_ratio() * g.intensity
            for g in source_lp_result_0.max_log_likelihood_fit
            .tracer_linear_light_profiles_to_light_profiles
            .galaxies[i + 1].bulge.profile_list
        ]
        total_luminosity = np.sum(luminosity_per_gaussian_list) / pixel_scale ** 2
        mass.einstein_radius = af.UniformPrior(
            lower_limit=0.0,
            upper_limit=min(5 * 0.5 * total_luminosity ** 0.6, 5.0),
        )

        extra_galaxies_fixed_list.append(af.Model(
            al.Galaxy, redshift=redshift_lens, bulge=extra_galaxy_bulge, mass=mass
        ))

    for i in range(len(outer_extra_galaxies_centres)):
        extra_galaxy_bulge = source_lp_result_0.instance.extra_galaxies[
            i + len(inner_extra_galaxies_centres)
        ].bulge

        mass = af.Model(al.mp.Isothermal)
        mass.centre = extra_galaxy_bulge.centre
        mass.ell_comps = extra_galaxy_bulge.ell_comps

        luminosity_per_gaussian_list = [
            2 * np.pi * g.sigma ** 2 / g.axis_ratio() * g.intensity
            for g in source_lp_result_0.max_log_likelihood_fit
            .tracer_linear_light_profiles_to_light_profiles
            .galaxies[i + 1 + len(inner_extra_galaxies_centres)].bulge.profile_list
        ]
        total_luminosity = np.sum(luminosity_per_gaussian_list) / pixel_scale ** 2
        mass.einstein_radius = scaling_factor * total_luminosity ** scaling_relation

        extra_galaxies_fixed_list.append(af.Model(
            al.Galaxy, redshift=redshift_lens, bulge=extra_galaxy_bulge, mass=mass
        ))

    extra_galaxies_fixed_centres = af.Collection(extra_galaxies_fixed_list)

    source_lp_result_1 = slam_pipeline.source_lp.run(
        settings_search=settings_search,
        analysis=analysis,
        lens_bulge=source_lp_result_0.instance.galaxies.lens.bulge,
        lens_disk=None,
        mass=lens_mass,
        shear=af.Model(al.mp.ExternalShear),
        source_bulge=source_bulge,
        extra_galaxies=extra_galaxies_fixed_centres,
        mass_centre=source_lp_result_0.instance.galaxies.lens.bulge.centre,
        redshift_lens=redshift_lens,
        redshift_source=redshift_source,
    )

    """
    __SOURCE PIX PIPELINE 1__

    The first SOURCE PIX search initializes the pixelized source reconstruction using
    a Delaunay mesh with AdaptSplit regularization. The adapt image is built from the
    source image estimated by the SOURCE LP PIPELINE.
    """
    hilbert_pixels = al.model_util.hilbert_pixels_from_pixel_scale(pixel_scale)

    image_mesh = al.image_mesh.Hilbert(pixels=hilbert_pixels, weight_power=3.5, weight_floor=0.01)

    galaxy_image_name_dict = al.galaxy_name_image_dict_via_result_from(
        result=source_lp_result_1
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
            source_lp_result_1.positions_likelihood_from(
                factor=2.0, positions=positions, minimum_threshold=0.2
            )
        ],
    )

    source_pix_result_1 = slam_pipeline.source_pix.run_1(
        settings_search=settings_search,
        analysis=analysis,
        source_lp_result=source_lp_result_1,
        extra_galaxies=extra_galaxies_fixed_centres,
        mesh_init=al.mesh.Delaunay(pixels=hilbert_pixels, zeroed_pixels=edge_pixels_total),
        regularization_init=af.Model(al.reg.AdaptSplit),
    )

    """
    __SOURCE PIX PIPELINE 2__

    The second SOURCE PIX search uses a Hilbert image mesh adapted to the source
    morphology and refines the mesh and regularization.
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
    )

    source_pix_result_2 = slam_pipeline.source_pix.run_2(
        settings_search=settings_search,
        analysis=analysis,
        source_lp_result=source_lp_result_1,
        source_pix_result_1=source_pix_result_1,
        mesh=al.mesh.Delaunay(pixels=hilbert_pixels, zeroed_pixels=edge_pixels_total),
        regularization=af.Model(al.reg.AdaptSplit),
    )

    """
    __LIGHT LP PIPELINE__

    The lens light is modeled with two independent MGE bases whose centres are free
    within ±0.2 arcsec of the lens centre. Using independent centre priors per basis
    (rather than a shared prior) allows the two components to describe spatially
    distinct light structures, which is important for group-scale lenses.

    Inner extra galaxy light profiles are re-fit with free centres; their mass is fixed
    from source_pix_result_1. Outer extra galaxy light and mass are both fixed.
    """
    centre_0_basis = [
        af.UniformPrior(lower_limit=lens_centre[0] - 0.2, upper_limit=lens_centre[0] + 0.2),
        af.UniformPrior(lower_limit=lens_centre[0] - 0.2, upper_limit=lens_centre[0] + 0.2),
    ]
    centre_1_basis = [
        af.UniformPrior(lower_limit=lens_centre[1] - 0.2, upper_limit=lens_centre[1] + 0.2),
        af.UniformPrior(lower_limit=lens_centre[1] - 0.2, upper_limit=lens_centre[1] + 0.2),
    ]

    total_gaussians = 30
    log10_sigma_list = np.linspace(-2, np.log10(mask_radius), total_gaussians)
    bulge_gaussian_list = []

    for j in range(2):
        gaussian_list = af.Collection(
            af.Model(al.lp_linear.Gaussian) for _ in range(total_gaussians)
        )
        for i, gaussian in enumerate(gaussian_list):
            gaussian.centre.centre_0 = centre_0_basis[j]
            gaussian.centre.centre_1 = centre_1_basis[j]
            gaussian.ell_comps = gaussian_list[0].ell_comps
            gaussian.sigma = 10 ** log10_sigma_list[i]
        bulge_gaussian_list += gaussian_list

    lens_bulge = af.Model(al.lp_basis.Basis, profile_list=bulge_gaussian_list)

    extra_galaxies_list = []

    for i, extra_galaxy_centre in enumerate(inner_extra_galaxies_centres):
        extra_galaxy_bulge = al.model_util.mge_model_from(
            mask_radius=mask_radius,
            total_gaussians=10,
            centre_prior_is_uniform=False,
            centre=(extra_galaxy_centre[0], extra_galaxy_centre[1]),
            centre_sigma=0.1,
        )
        mass = source_pix_result_1.instance.extra_galaxies[i].mass
        extra_galaxies_list.append(af.Model(
            al.Galaxy, redshift=redshift_lens, bulge=extra_galaxy_bulge, mass=mass
        ))

    for i in range(len(outer_extra_galaxies_centres)):
        bulge = source_pix_result_1.instance.extra_galaxies[
            i + len(inner_extra_galaxies_centres)
        ].bulge
        mass = source_pix_result_1.instance.extra_galaxies[
            i + len(inner_extra_galaxies_centres)
        ].mass
        extra_galaxies_list.append(af.Model(
            al.Galaxy, redshift=redshift_lens, bulge=bulge, mass=mass
        ))

    extra_galaxies = af.Collection(extra_galaxies_list)

    light_result_1 = slam_pipeline.light_lp.run(
        settings_search=settings_search,
        analysis=analysis,
        source_result_for_lens=source_pix_result_1,
        source_result_for_source=source_pix_result_2,
        extra_galaxies=extra_galaxies,
        lens_bulge=lens_bulge,
        lens_disk=None,
    )

    """
    __MASS TOTAL PIPELINE__

    Inner extra galaxy mass profiles are re-fit with Einstein radii bounded by twice
    their light_result_1 values. Outer extra galaxy masses use a fresh luminosity
    scaling relation computed from the light_result_1 fit.
    """
    analysis = al.AnalysisImaging(
        dataset=dataset,
        adapt_images=adapt_images,
        positions_likelihood_list=[
            source_lp_result_1.positions_likelihood_from(
                factor=2.0, positions=positions, minimum_threshold=0.2
            )
        ],
    )

    extra_galaxies_list = []

    for i in range(len(inner_extra_galaxies_centres)):
        extra_galaxy_bulge = light_result_1.instance.extra_galaxies[i].bulge

        mass = af.Model(al.mp.Isothermal)
        mass.centre = extra_galaxy_bulge.centre
        mass.ell_comps = extra_galaxy_bulge.ell_comps
        mass.einstein_radius = af.UniformPrior(
            lower_limit=0.0,
            upper_limit=2 * light_result_1.instance.extra_galaxies[i].mass.einstein_radius,
        )

        extra_galaxies_list.append(af.Model(
            al.Galaxy, redshift=redshift_lens, bulge=extra_galaxy_bulge, mass=mass
        ))

    scaling_factor = af.UniformPrior(lower_limit=0.0, upper_limit=0.5)
    scaling_relation = af.UniformPrior(lower_limit=0.0, upper_limit=2.0)

    for i in range(len(outer_extra_galaxies_centres)):
        extra_galaxy_bulge = light_result_1.instance.extra_galaxies[
            i + len(inner_extra_galaxies_centres)
        ].bulge

        mass = af.Model(al.mp.Isothermal)
        mass.centre = extra_galaxy_bulge.centre
        mass.ell_comps = extra_galaxy_bulge.ell_comps

        luminosity_per_gaussian_list = [
            2 * np.pi * g.sigma ** 2 / g.axis_ratio() * g.intensity
            for g in light_result_1.max_log_likelihood_fit
            .tracer_linear_light_profiles_to_light_profiles
            .galaxies[i + 1 + len(inner_extra_galaxies_centres)].bulge.profile_list
        ]
        total_luminosity = np.sum(luminosity_per_gaussian_list) / pixel_scale ** 2
        mass.einstein_radius = scaling_factor * total_luminosity ** scaling_relation

        extra_galaxies_list.append(af.Model(
            al.Galaxy, redshift=redshift_lens, bulge=extra_galaxy_bulge, mass=mass
        ))

    extra_galaxies = af.Collection(extra_galaxies_list)

    mass_result = slam_pipeline.mass_total.run(
        settings_search=settings_search,
        analysis=analysis,
        source_result_for_lens=source_pix_result_1,
        source_result_for_source=source_pix_result_2,
        light_result=light_result_1,
        mass=af.Model(al.mp.PowerLaw),
        extra_galaxies=extra_galaxies,
        multipole_1=af.Model(al.mp.PowerLawMultipole),
        multipole_3=af.Model(al.mp.PowerLawMultipole),
        multipole_4=af.Model(al.mp.PowerLawMultipole),
    )

    return source_lp_result_1, mass_result


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="PyAutoLens SLAM Group Pipeline")

    parser.add_argument(
        "--sample",
        metavar="name",
        required=False,
        default=None,
        help="Name of the sample subdirectory inside dataset/ (e.g. euclid_groups).",
    )

    parser.add_argument(
        "--dataset",
        metavar="name",
        required=False,
        default=None,
        help="Name of the dataset subdirectory inside dataset/[sample]/ (e.g. group0001).",
    )

    args = parser.parse_args()

    fit(dataset_name=args.dataset, sample_name=args.sample)


"""
Finish.
"""
