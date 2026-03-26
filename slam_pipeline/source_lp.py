import autofit as af
import autolens as al

from typing import Union, Optional, Tuple


def run(
    settings_search: af.SettingsSearch,
    analysis: Union[al.AnalysisImaging, al.AnalysisInterferometer],
    lens_bulge: Optional[af.Model] = af.Model(al.lp.Sersic),
    lens_disk: Optional[af.Model] = af.Model(al.lp.Exponential),
    lens_point: Optional[af.Model] = None,
    mass: af.Model = af.Model(al.mp.Isothermal),
    shear: af.Model(al.mp.ExternalShear) = af.Model(al.mp.ExternalShear),
    source_bulge: Optional[af.Model] = af.Model(al.lp.Sersic),
    source_disk: Optional[af.Model] = None,
    redshift_lens: float = 0.5,
    redshift_source: float = 1.0,
    mass_centre: Optional[Tuple[float, float]] = None,
    extra_galaxies: Optional[af.Collection] = None,
    dataset_model: Optional[af.Model] = None,
) -> af.Result:
    """
    The SlaM SOURCE LP PIPELINE, which provides an initial model for the lens's light, mass and source using a
    parametric source model (e.g. MGE, Sersics).

    Parameters
    ----------
    analysis
        The analysis class which includes the `log_likelihood_function` and can be customized for the SLaM model-fit.
    lens_bulge
        The model used to represent the light distribution of the lens galaxy's bulge (set to
        None to omit a bulge).
    lens_disk
        The model used to represent the light distribution of the lens galaxy's disk (set to
        None to omit a disk).
    lens_point
        The model used to represent the light distribution of point source emission in the lens galaxy (set to
        None to omit a point).
    mass
        The `MassProfile` fitted by this pipeline.
    shear
        The model used to represent the external shear in the mass model (set to None to turn off shear).
    source_bulge
        The model used to represent the light distribution of the source galaxy's bulge (set to
        None to omit a bulge).
    source_disk
        The model used to represent the light distribution of the source galaxy's disk (set to
        None to omit a disk).
    redshift_lens
        The redshift of the lens galaxy fitted, used by the pipeline for converting arc-seconds to kpc, masses to
        solMass, etc.
    redshift_source
        The redshift of the source galaxy fitted, used by the pipeline for converting arc-seconds to kpc, masses to
        solMass, etc.
    mass_centre
       If input, a fixed (y,x) centre of the mass profile is used which is not treated as a free parameter by the
       non-linear search.
    extra_galaxies
        Additional extra galaxies containing light and mass profiles, which model nearby line of sight galaxies.
    dataset_model
        Add aspects of the dataset to the model, for example the arc-second (y,x) offset between two datasets for
        multi-band fitting or the background sky level.
    """

    """
    __Model + Search + Analysis + Model-Fit (Search 1)__

    Search 1 of the SOURCE LP PIPELINE fits a lens model where:

     - The lens galaxy light is modeled using a light profiles [no prior initialization].
     - The lens galaxy mass is modeled using a total mass distribution [no prior initialization].
     - The source galaxy's light is a light profiles [no prior initialization].

    This search aims to accurately estimate an initial lens light model, mass model and source model.
    """

    if mass_centre is not None:
        mass.centre = mass_centre

    model = af.Collection(
        galaxies=af.Collection(
            lens=af.Model(
                al.Galaxy,
                redshift=redshift_lens,
                bulge=lens_bulge,
                disk=lens_disk,
                point=lens_point,
                mass=mass,
                shear=shear,
            ),
            source=af.Model(
                al.Galaxy,
                redshift=redshift_source,
                bulge=source_bulge,
                disk=source_disk,
            ),
        ),
        extra_galaxies=extra_galaxies,
        dataset_model=dataset_model,
    )

    search = af.Nautilus(
        name="source_lp[1]",
        **settings_search.search_dict,
        n_live=150,
        n_batch=50,
        n_like_max=200000
    )

    result = search.fit(model=model, analysis=analysis, **settings_search.fit_dict)

    return result


def run_0__group(
    settings_search: af.SettingsSearch,
    analysis: Union[al.AnalysisImaging, al.AnalysisInterferometer],
    lens_galaxy_models: list,
    extra_galaxies: Optional[af.Collection] = None,
    scaling_galaxies: Optional[af.Collection] = None,
    dataset_model: Optional[af.Model] = None,
) -> af.Result:
    """
    An initial SOURCE LP search used for group-scale lenses, which fits only the lens
    light profiles (one per main lens galaxy) with no mass model or source galaxy.

    This provides a clean initialization of all light components — and measured
    luminosities for each extra and scaling galaxy — before the mass model and source
    are introduced in `run_group`.

    Parameters
    ----------
    settings_search
        AutoFit settings controlling output paths and search configuration.
    analysis
        The analysis class containing the log-likelihood function.
    lens_galaxy_models
        A list of ``af.Model(al.Galaxy)`` instances, one per main lens galaxy, each
        containing only light profiles (bulge/disk/point) with no mass component.
        Built by the caller so that centres and priors can be set per-lens.
    extra_galaxies
        Nearby extra galaxies with free-centre light profiles.
    scaling_galaxies
        More distant galaxies with free-centre light profiles whose luminosities are
        later used to set Einstein radii via a shared scaling relation.  Kept separate
        from ``extra_galaxies`` so the two populations remain distinguishable when
        building stage-1 models.
    """
    lens_dict = {f"lens_{i}": m for i, m in enumerate(lens_galaxy_models)}

    model = af.Collection(
        galaxies=af.Collection(**lens_dict),
        extra_galaxies=extra_galaxies,
        scaling_galaxies=scaling_galaxies,
        dataset_model=dataset_model,
    )

    n_extra = len(extra_galaxies) if extra_galaxies is not None else 0
    n_scaling = len(scaling_galaxies) if scaling_galaxies is not None else 0
    n_live = 100 + 30 * len(lens_galaxy_models) + 30 * n_extra + 30 * n_scaling

    search = af.Nautilus(
        name="source_lp[0]",
        **settings_search.search_dict,
        n_live=n_live,
        n_batch=50,
        n_like_max=1000000,
    )

    return search.fit(model=model, analysis=analysis, **settings_search.fit_dict)


def run_group(
    settings_search: af.SettingsSearch,
    analysis: Union[al.AnalysisImaging, al.AnalysisInterferometer],
    lens_galaxy_models: list,
    source_bulge: Optional[af.Model] = af.Model(al.lp.Sersic),
    source_disk: Optional[af.Model] = None,
    redshift_source: float = 1.0,
    extra_galaxies: Optional[af.Collection] = None,
    scaling_galaxies: Optional[af.Collection] = None,
    dataset_model: Optional[af.Model] = None,
) -> af.Result:
    """
    The SLaM SOURCE LP PIPELINE for group-scale lenses with multiple main lens galaxies.

    Fits all main lens light and mass profiles simultaneously with a parametric source.
    This is the group-scale counterpart to ``run``, where ``galaxies`` contains
    ``lens_0``, ``lens_1``, ... rather than a single ``lens``.

    Typically called after ``run_0__group`` so that the light profiles in
    ``lens_galaxy_models`` can be initialized from that result by the caller.

    Parameters
    ----------
    settings_search
        AutoFit settings controlling output paths and search configuration.
    analysis
        The analysis class containing the log-likelihood function.
    lens_galaxy_models
        A list of ``af.Model(al.Galaxy)`` instances, one per main lens galaxy,
        each containing bulge/disk/point, mass, and shear components as appropriate.
        Built and configured by the caller (centres, priors, prior initialization
        from ``run_0__group`` result).
    source_bulge
        The model for the source galaxy's bulge light.
    source_disk
        The model for the source galaxy's disk light (``None`` to omit).
    redshift_source
        The redshift of the source galaxy.
    extra_galaxies
        Nearby extra galaxies modeled with individual light and mass profiles
        (renamed from ``inner_extra_galaxies``).
    scaling_galaxies
        More distant extra galaxies modeled with a shared luminosity scaling
        relation for their mass (renamed from ``outer_extra_galaxies``).
    dataset_model
        Optional dataset-level model components (e.g. background sky, astrometric
        offsets for multi-band fitting).
    """
    lens_dict = {f"lens_{i}": m for i, m in enumerate(lens_galaxy_models)}
    lens_dict["source"] = af.Model(
        al.Galaxy,
        redshift=redshift_source,
        bulge=source_bulge,
        disk=source_disk,
    )


    model = af.Collection(
        galaxies=af.Collection(**lens_dict),
        extra_galaxies=extra_galaxies,
        scaling_galaxies=scaling_galaxies,
        dataset_model=dataset_model,
    )

    n_extra = len(extra_galaxies) if extra_galaxies is not None else 0
    n_scaling = len(scaling_galaxies) if scaling_galaxies is not None else 0
    n_live = 150 + 30 * len(lens_galaxy_models) + 30 * n_extra + 30 * n_scaling

    search = af.Nautilus(
        name="source_lp[1]",
        **settings_search.search_dict,
        n_live=n_live,
        n_batch=50,
        n_like_max=200000,
    )

    return search.fit(model=model, analysis=analysis, **settings_search.fit_dict)
