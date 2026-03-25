"""
SLaM (Source, Light and Mass): Group — Multiple Main Lenses
============================================================

This script extends the group pipeline to support multiple main lens galaxies
(e.g. a merging pair or a group-core triplet). The number of main lenses is
determined at runtime from ``main_lens_centres.json``.

Galaxy populations
------------------
- **Main lenses** (``lens_0``, ``lens_1``, ...): each has an MGE bulge, an
  Isothermal mass, and ExternalShear. Centres are pinned to the positions in
  ``main_lens_centres.json``.
- **Extra lens galaxies**: nearby companions with individual MGE light + Isothermal
  mass. Loaded from ``extra_galaxies_centres.json`` (optional).
- **Scaling lens galaxies**: more distant galaxies modeled with a shared
  luminosity scaling relation. Loaded from ``scaling_galaxies_centres.json``
  (optional).

Pipeline stages implemented here
---------------------------------
source_lp[0]  run_0__group  — light only, no mass, no source
source_lp[1]  run_group     — full LP fit, mass + source introduced

Subsequent stages (source_pix, light_lp, mass_total) will be added once the
group variants of those pipeline functions are implemented.
"""
import numpy as np
import sys
import json
from pathlib import Path


def _load_centres(path):
    """Load a centres JSON file, returning an empty list if the file is absent."""
    try:
        import autolens as al
        return al.Grid2DIrregular(al.from_json(file_path=path))
    except FileNotFoundError:
        import autolens as al
        return al.Grid2DIrregular([])


def fit(dataset_name, sample_name):
    """
    __Paths__
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
    __Galaxy Centres__

    main_lens_centres.json        — required; determines the number of main lenses.
    extra_galaxies_centres.json  — optional; empty list if absent.
    scaling_galaxies_centres.json — optional; empty list if absent.

    All three files contain a list of [y, x] arcsecond coordinates.
    """
    main_lens_centres = _load_centres(dataset_path / "main_lens_centres.json")
    extra_lens_centres = _load_centres(dataset_path / "extra_galaxies_centres.json")
    scaling_lens_centres = _load_centres(dataset_path / "scaling_galaxies_centres.json")

    all_galaxy_centres = al.Grid2DIrregular(
        main_lens_centres.in_list
        + extra_lens_centres.in_list
        + scaling_lens_centres.in_list
    )

    positions = al.Grid2DIrregular(
        al.from_json(file_path=dataset_path / "positions.json")
    )

    """
    __Mask__
    """
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
        centre_list=list(all_galaxy_centres),
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
    __SOURCE LP PIPELINE — stage 0: light only__

    Fits MGE light profiles for every main lens galaxy simultaneously with no
    mass model and no source. Extra lens galaxies are also fit with light-only
    MGE profiles.

    The result initializes:
      - bulge centres for each main lens (used as mass_centre in stage 1)
      - extra lens galaxy luminosities (used to bound mass Einstein radii in stage 1)
    """
    analysis = al.AnalysisImaging(dataset=dataset)

    # --- main lens light models (one per centre, light only) ---
    lens_light_models = []
    for centre in main_lens_centres:
        bulge = al.model_util.mge_model_from(
            mask_radius=mask_radius,
            total_gaussians=30,
            gaussian_per_basis=2,
            centre_prior_is_uniform=False,
            centre=(centre[0], centre[1]),
            centre_sigma=0.1,
        )
        lens_light_models.append(
            af.Model(al.Galaxy, redshift=redshift_lens, bulge=bulge)
        )

    # --- extra lens galaxy light models ---
    extra_lens_light_list = []
    for centre in extra_lens_centres:
        bulge = al.model_util.mge_model_from(
            mask_radius=mask_radius,
            total_gaussians=10,
            centre_prior_is_uniform=False,
            centre=(centre[0], centre[1]),
            centre_sigma=0.1,
        )
        extra_lens_light_list.append(
            af.Model(al.Galaxy, redshift=redshift_lens, bulge=bulge)
        )

    extra_galaxies_free = (
        af.Collection(extra_lens_light_list) if extra_lens_light_list else None
    )

    source_lp_result_0 = slam_pipeline.source_lp.run_0__group(
        settings_search=settings_search,
        analysis=analysis,
        lens_galaxy_models=lens_light_models,
        extra_galaxies=extra_galaxies_free,
    )

    """
    __SOURCE LP PIPELINE — stage 1: full LP fit__

    Introduces mass and source. For each main lens:
      - bulge is fixed to the stage-0 instance (light already well-characterized)
      - mass centre is pinned to the stage-0 bulge centre
      - mass Einstein radius has a broad uniform prior
      - ExternalShear is free

    Extra lens galaxies have light fixed from stage 0; mass Einstein radii are
    bounded by the stage-0 luminosity (same luminosity-scaling bound as group.py).

    Scaling lens galaxies have mass modeled via a shared luminosity scaling relation.
    """
    analysis = al.AnalysisImaging(
        dataset=dataset,
        positions_likelihood_list=[
            al.PositionsLH(positions=positions, threshold=0.3)
        ],
    )

    source_bulge = al.model_util.mge_model_from(
        mask_radius=1.0,
        total_gaussians=30,
        centre_prior_is_uniform=False,
        centre=(main_lens_centres[0][0], main_lens_centres[0][1]),
        centre_sigma=0.6,
    )

    # --- main lens full models (light fixed from stage 0, mass + shear free) ---
    lens_full_models = []
    for i in range(len(main_lens_centres)):
        lp0_lens = getattr(source_lp_result_0.instance.galaxies, f"lens_{i}")

        mass = af.Model(al.mp.Isothermal)
        mass.centre = lp0_lens.bulge.centre
        mass.einstein_radius = af.UniformPrior(lower_limit=0.0, upper_limit=5.0)

        lens_full_models.append(
            af.Model(
                al.Galaxy,
                redshift=redshift_lens,
                bulge=lp0_lens.bulge,
                mass=mass,
                shear=af.Model(al.mp.ExternalShear),
            )
        )

    # --- extra lens galaxy models (light fixed, mass bounded by luminosity) ---
    extra_lens_fixed_list = []
    for i in range(len(extra_lens_centres)):
        lp0_extra = source_lp_result_0.instance.extra_galaxies[i]

        mass = af.Model(al.mp.Isothermal)
        mass.centre = lp0_extra.bulge.centre
        mass.ell_comps = lp0_extra.bulge.ell_comps

        luminosity_per_gaussian_list = [
            2 * np.pi * g.sigma ** 2 / g.axis_ratio() * g.intensity
            for g in (
                source_lp_result_0.max_log_likelihood_fit
                .tracer_linear_light_profiles_to_light_profiles
                .extra_galaxies[i].bulge.profile_list
            )
        ]
        total_luminosity = np.sum(luminosity_per_gaussian_list) / pixel_scale ** 2
        mass.einstein_radius = af.UniformPrior(
            lower_limit=0.0,
            upper_limit=min(5 * 0.5 * total_luminosity ** 0.6, 5.0),
        )

        extra_lens_fixed_list.append(
            af.Model(al.Galaxy, redshift=redshift_lens, bulge=lp0_extra.bulge, mass=mass)
        )

    extra_galaxies_fixed = (
        af.Collection(extra_lens_fixed_list) if extra_lens_fixed_list else None
    )

    # --- scaling lens galaxy models (shared luminosity scaling relation) ---
    scaling_factor = af.UniformPrior(lower_limit=0.0, upper_limit=0.5)
    scaling_relation = af.UniformPrior(lower_limit=0.0, upper_limit=2.0)

    scaling_lens_list = []
    for i in range(len(scaling_lens_centres)):
        # Scaling galaxies had no light in stage 0; use their known centre directly.
        centre = scaling_lens_centres[i]

        mass = af.Model(al.mp.Isothermal)
        mass.centre = (centre[0], centre[1])
        # Einstein radius driven by a shared scaling relation; no luminosity bound.
        mass.einstein_radius = scaling_factor * scaling_relation

        scaling_lens_list.append(
            af.Model(al.Galaxy, redshift=redshift_lens, mass=mass)
        )

    scaling_galaxies = (
        af.Collection(scaling_lens_list) if scaling_lens_list else None
    )

    source_lp_result_1 = slam_pipeline.source_lp.run_group(
        settings_search=settings_search,
        analysis=analysis,
        lens_galaxy_models=lens_full_models,
        source_bulge=source_bulge,
        redshift_source=redshift_source,
        extra_galaxies=extra_galaxies_fixed,
        scaling_galaxies=scaling_galaxies,
    )

    """
    __SUBSEQUENT STAGES__

    source_pix.run_1_group, source_pix.run_2_group, light_lp.run_group, and
    mass_total.run_group will be added here once those pipeline functions are
    implemented.

    Result access pattern for downstream stages:

        n_lenses = sum(
            1 for k in vars(source_lp_result_1.instance.galaxies)
            if k.startswith("lens_")
        )
        for i in range(n_lenses):
            lp1_lens = getattr(source_lp_result_1.instance.galaxies, f"lens_{i}")
            lp1_lens_model = getattr(source_lp_result_1.model.galaxies, f"lens_{i}")
    """

    return source_lp_result_0, source_lp_result_1


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="PyAutoLens SLAM Group Pipeline — Multiple Main Lenses"
    )

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
