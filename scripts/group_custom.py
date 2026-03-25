"""
SLaM (Source, Light and Mass): Group — Custom
==============================================

Extends `group.py` with a source-only image extraction step performed immediately
after `source_lp[1]` completes.

__Source Only Image__

After the SOURCE LP PIPELINE has a clean model for the lensed source, a
source-only version of the dataset is written to disk. The procedure is:

1. Extract the S/N map of the model source image from the source_lp[1] result.
2. Build a binary mask of pixels where source S/N > 1.5.
3. For pixels **outside** the source mask, replace the image value with a random
   Gaussian draw N(0, σ) where σ is the actual per-pixel noise, so the replaced
   pixels look statistically like blank sky.
4. Construct a new ``al.Imaging`` object directly from the modified arrays, apply
   ``apply_noise_scaling`` to inflate the noise at non-source pixels, then apply
   the same circular mask and over-sampling as the main dataset.

Everything else is identical to ``group.py``.
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
    __Hilbert Pixels__
    """
    hilbert_pixels = al.model_util.hilbert_pixels_from_pixel_scale(pixel_scale)

    """
    __SOURCE LP PIPELINE__

    The SOURCE LP PIPELINE runs in two stages.

    The first search (run_0__group) fits only the lens and extra galaxy light profiles
    with no mass model. This initializes the light model and measures the extra galaxy
    luminosities before the mass model is introduced.

    The second search (run) introduces the lens mass and source light, holding the extra
    galaxy light profiles fixed to their run_0__group values. Inner extra galaxy mass
    profiles are initialized with Einstein radii bounded by their run_0__group
    luminosities; outer extra galaxy mass profiles use a shared luminosity scaling
    relation.
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
    __Source Only Image__

    The source_lp[1] result provides a clean model of the lensed source. We use
    this to build a source-only version of the dataset in memory:

    1. Extract the S/N map for the source galaxy from the source_lp[1] result.
       Pixels where S/N > 1.5 are flagged as containing source signal.

    2. Pixels **outside** the source signal mask are replaced with a random
       Gaussian draw N(0, σ), where σ is the actual per-pixel noise, so the
       replaced sky looks statistically identical to blank data.

    3. A new ``al.Imaging`` object is constructed from the modified arrays
       (reusing the original PSF). ``apply_noise_scaling`` inflates the noise at
       non-source pixels so they contribute negligibly to the likelihood.
    """
    galaxy_image_name_dict = al.galaxy_name_image_dict_via_result_from(
        result=source_lp_result_1
    )

    # galaxy_name_image_dict_via_result_from returns S/N maps (not flux), so the
    # values are already signal-to-noise — no division by noise_map is needed.
    source_snr = galaxy_image_name_dict["('galaxies', 'source')"]
    source_signal_mask = source_snr.native > 1.5

    data_native = dataset.data.native
    noise_native = dataset.noise_map.native

    # Draw blank-sky noise for pixels outside the source mask. Use the local
    # per-pixel noise level, falling back to 1.0 for any zero-valued pixels
    # (which lie outside the circular mask and will be re-masked on reload).
    noise_scale = np.where(noise_native > 0, noise_native, 1.0)
    rng = np.random.default_rng(seed=0)
    data_source_only = np.where(
        source_signal_mask,
        data_native,
        rng.normal(loc=0.0, scale=noise_scale),
    )

    dataset_source_only = al.Imaging(
        data=al.Array2D.no_mask(values=data_source_only, pixel_scales=pixel_scale),
        noise_map=al.Array2D.no_mask(values=noise_scale, pixel_scales=pixel_scale),
        psf=dataset.psf,
    )

    source_outside_mask = al.Mask2D(
        mask=~source_signal_mask, pixel_scales=pixel_scale
    )
    dataset_source_only = dataset_source_only.apply_noise_scaling(mask=source_outside_mask)
    dataset_source_only = dataset_source_only.apply_mask(mask=mask)
    dataset_source_only = dataset_source_only.apply_over_sampling(
        over_sample_size_lp=over_sample_size
    )

    settings_search_source_only = af.SettingsSearch(
        path_prefix=sample_name,
        unique_tag=f"{dataset_name}_source_only",
        info=info,
        session=None,
    )


    """
    __SOURCE-ONLY SOURCE PIX PIPELINE 1__

    Run the first pixelized source pipeline on the source-only dataset with no lens
    light model (lens_bulge=None). Mass and shear are initialized from
    source_lp_result_1. Adapt images are the S/N maps already computed from
    source_lp_result_1.
    """
    edge_pixels_total = 30
    signal_to_noise_threshold = 3.0

    image_mesh_so = al.image_mesh.Hilbert(pixels=hilbert_pixels, weight_power=3.5, weight_floor=0.01)

    image_plane_mesh_grid_so = image_mesh_so.image_plane_mesh_grid_from(
        mask=dataset_source_only.mask,
        adapt_data=galaxy_image_name_dict["('galaxies', 'source')"],
    )

    image_plane_mesh_grid_so = al.image_mesh.append_with_circle_edge_points(
        image_plane_mesh_grid=image_plane_mesh_grid_so,
        centre=mask.mask_centre,
        radius=mask_radius + mask.pixel_scale / 2.0,
        n_points=edge_pixels_total,
    )

    adapt_images_so = al.AdaptImages(
        galaxy_name_image_dict=galaxy_image_name_dict,
        galaxy_name_image_plane_mesh_grid_dict={
            "('galaxies', 'source')": image_plane_mesh_grid_so
        },
    )

    over_sample_size_pixelization_so = np.where(
        galaxy_image_name_dict["('galaxies', 'source')"] > signal_to_noise_threshold,
        4,
        2,
    )
    over_sample_size_pixelization_so = al.Array2D(
        values=over_sample_size_pixelization_so, mask=mask
    )

    dataset_source_only = dataset_source_only.apply_over_sampling(
        over_sample_size_lp=over_sample_size,
        over_sample_size_pixelization=over_sample_size_pixelization_so,
    )

    mass_so = al.util.chaining.mass_from(
        mass=source_lp_result_1.model.galaxies.lens.mass,
        mass_result=source_lp_result_1.model.galaxies.lens.mass,
        unfix_mass_centre=True,
    )

    analysis_source_only = al.AnalysisImaging(
        dataset=dataset_source_only,
        adapt_images=adapt_images_so,
        positions_likelihood_list=[
            source_lp_result_1.positions_likelihood_from(
                factor=2.0, positions=positions, minimum_threshold=0.2
            )
        ],
    )

    source_pix_so_result_1 = slam_pipeline.source_pix.run_1__bypass_lp(
        settings_search=settings_search_source_only,
        analysis=analysis_source_only,
        lens_bulge=None,
        lens_disk=None,
        mass=mass_so,
        shear=source_lp_result_1.model.galaxies.lens.shear,
        extra_galaxies=extra_galaxies_fixed_centres,
        mesh_init=al.mesh.Delaunay(pixels=hilbert_pixels, zeroed_pixels=edge_pixels_total),
        regularization_init=af.Model(al.reg.AdaptSplit),
        redshift_lens=redshift_lens,
        redshift_source=redshift_source,
    )

    """
    __SOURCE-ONLY SOURCE PIX PIPELINE 2__

    Refine the pixelized source with adapt images from source_pix_so_result_1.
    """
    galaxy_image_name_dict_so = al.galaxy_name_image_dict_via_result_from(
        result=source_pix_so_result_1
    )

    image_mesh_so = al.image_mesh.Hilbert(pixels=hilbert_pixels, weight_power=3.5, weight_floor=0.01)

    signal_to_noise_threshold_so = 3.0
    adapt_data_snr_max_so = galaxy_image_name_dict_so["('galaxies', 'source')"]
    adapt_data_snr_max_so[adapt_data_snr_max_so > signal_to_noise_threshold_so] = (
        signal_to_noise_threshold_so
    )

    image_plane_mesh_grid_so = image_mesh_so.image_plane_mesh_grid_from(
        mask=dataset_source_only.mask, adapt_data=adapt_data_snr_max_so
    )

    image_plane_mesh_grid_so = al.image_mesh.append_with_circle_edge_points(
        image_plane_mesh_grid=image_plane_mesh_grid_so,
        centre=mask.mask_centre,
        radius=mask_radius + mask.pixel_scale / 2.0,
        n_points=edge_pixels_total,
    )

    adapt_images_so = al.AdaptImages(
        galaxy_name_image_dict=galaxy_image_name_dict_so,
        galaxy_name_image_plane_mesh_grid_dict={
            "('galaxies', 'source')": image_plane_mesh_grid_so
        },
    )

    over_sample_size_pixelization_so = np.where(
        galaxy_image_name_dict_so["('galaxies', 'source')"] > signal_to_noise_threshold,
        4,
        2,
    )
    over_sample_size_pixelization_so = al.Array2D(
        values=over_sample_size_pixelization_so, mask=mask
    )

    dataset_source_only = dataset_source_only.apply_over_sampling(
        over_sample_size_lp=over_sample_size,
        over_sample_size_pixelization=over_sample_size_pixelization_so,
    )

    analysis_source_only = al.AnalysisImaging(
        dataset=dataset_source_only,
        adapt_images=adapt_images_so,
    )

    source_pix_so_result_2 = slam_pipeline.source_pix.run_2(
        settings_search=settings_search_source_only,
        analysis=analysis_source_only,
        source_lp_result=source_pix_so_result_1,
        source_pix_result_1=source_pix_so_result_1,
        mesh=al.mesh.Delaunay(pixels=hilbert_pixels, zeroed_pixels=edge_pixels_total),
        regularization=af.Model(al.reg.AdaptSplit),
    )

    """
    __SOURCE-ONLY MASS TOTAL PIPELINE__

    Fit the mass model on the source-only dataset with no lens light (light_result=None).
    Extra galaxy masses are fixed to the values from source_pix_so_result_1.
    """
    analysis_source_only = al.AnalysisImaging(
        dataset=dataset_source_only,
        adapt_images=adapt_images_so,
        positions_likelihood_list=[
            source_lp_result_1.positions_likelihood_from(
                factor=2.0, positions=positions, minimum_threshold=0.2
            )
        ],
    )

    # Extra Galaxies:

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
            al.Galaxy, redshift=redshift_lens, mass=mass
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
            al.Galaxy, redshift=redshift_lens, mass=mass
        ))

    extra_galaxies_fixed_centres = af.Collection(extra_galaxies_fixed_list)

    mass_so_result = slam_pipeline.mass_total.run(
        settings_search=settings_search_source_only,
        analysis=analysis_source_only,
        source_result_for_lens=source_pix_so_result_1,
        source_result_for_source=source_pix_so_result_2,
        light_result=None,
        mass=af.Model(al.mp.PowerLaw),
        extra_galaxies=extra_galaxies_fixed_centres,
        multipole_1=af.Model(al.mp.PowerLawMultipole),
        multipole_3=af.Model(al.mp.PowerLawMultipole),
        multipole_4=af.Model(al.mp.PowerLawMultipole),
    )

    return source_lp_result_1, mass_so_result


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="PyAutoLens SLAM Group Pipeline (Custom)")

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
