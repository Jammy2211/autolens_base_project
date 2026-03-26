import autofit as af
import autolens as al

from typing import List, Union, Optional, Tuple


def run(
    settings_search: af.SettingsSearch,
    analysis: Union[al.AnalysisImaging, al.AnalysisInterferometer],
    source_result_for_lens: af.Result,
    source_result_for_source: af.Result,
    light_result: Optional[af.Result],
    mass: af.Model = af.Model(al.mp.Isothermal),
    multipole_1: Optional[af.Model] = None,
    multipole_3: Optional[af.Model] = None,
    multipole_4: Optional[af.Model] = None,
    smbh: Optional[af.Model] = None,
    mass_centre: Optional[Tuple[float, float]] = None,
    extra_galaxies: Optional[af.Collection] = None,
    dataset_model: Optional[af.Model] = None,
    n_batch: int = 20,
) -> af.Result:
    """
    The SLaM MASS TOTAL PIPELINE, which fits a lens model with a total mass distribution (e.g. a power-law).

    Parameters
    ----------
    analysis
        The analysis class which includes the `log_likelihood_function` and can be customized for the SLaM model-fit.
    source_result_for_lens
        The result of the SLaM SOURCE LP PIPELINE or SOURCE PIX PIPELINE which ran before this pipeline,
        used for initializing model components associated with the lens galaxy.
    source_result_for_source
        The result of the SLaM SOURCE LP PIPELINE or SOURCE PIX PIPELINE which ran before this pipeline,
        used for initializing model components associated with the source galaxy.
    light_result
        The result of the SLaM LIGHT LP PIPELINE which ran before this pipeline.
    mass
        The `MassProfile` used to fit the lens galaxy mass in this pipeline.
    light_linear_to_standard
        If `True`, convert all linear light profiles in the model to standard light profiles, whose `intensity` values
        use the max likelihood result of the LIGHT PIPELINE.
    multipole_1
        Optionally include a first order multipole mass profile component in the mass model.
    multipole_3
        Optionally include a third order multipole mass profile component in the mass model.
    multipole_4
        Optionally include a fourth order multipole mass profile component in the mass model.
    smbh
        The `MassProfile` used to fit the a super massive black hole in the lens galaxy.
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

    Search 1 of the MASS TOTAL PIPELINE fits a lens model where:

     - The lens galaxy mass is modeled using a total mass distribution [Priors initialized from SOURCE PIPELINE].
     - The source galaxy's light is parametric or a pixelization depending on the previous pipeline [Model and priors 
     initialized from SOURCE PIPELINE].

    This search aims to accurately estimate the lens mass model, using the improved mass model priors and source model 
    of the SOURCE PIPELINE
    """
    mass = al.util.chaining.mass_from(
        mass=mass,
        mass_result=source_result_for_lens.model.galaxies.lens.mass,
        unfix_mass_centre=True,
    )

    if mass_centre is not None:
        mass.centre = mass_centre

    if smbh is not None:
        smbh.centre = mass.centre

    if light_result is None:
        bulge = None
        disk = None
        point = None

    else:
        bulge = light_result.instance.galaxies.lens.bulge
        disk = light_result.instance.galaxies.lens.disk
        point = light_result.instance.galaxies.lens.point

    n_live_multipole = 0

    if multipole_1 is not None:
        multipole_1.m = 1
        multipole_1.centre = mass.centre
        multipole_1.einstein_radius = mass.einstein_radius
        multipole_1.slope = mass.slope

        n_live_multipole += 100

    if multipole_3 is not None:
        multipole_3.m = 3
        multipole_3.centre = mass.centre
        multipole_3.einstein_radius = mass.einstein_radius
        multipole_3.slope = mass.slope

        n_live_multipole += 100

    if multipole_4 is not None:
        multipole_4.m = 4
        multipole_4.centre = mass.centre
        multipole_4.einstein_radius = mass.einstein_radius
        multipole_4.slope = mass.slope

        n_live_multipole += 100

    source = al.util.chaining.source_from(
        result=source_result_for_source,
    )

    model = af.Collection(
        galaxies=af.Collection(
            lens=af.Model(
                al.Galaxy,
                redshift=source_result_for_lens.instance.galaxies.lens.redshift,
                bulge=bulge,
                disk=disk,
                point=point,
                mass=mass,
                multipole_1=multipole_1,
                multipole_3=multipole_3,
                multipole_4=multipole_4,
                shear=source_result_for_lens.model.galaxies.lens.shear,
                smbh=smbh,
            ),
            source=source,
        ),
        extra_galaxies=extra_galaxies,
        dataset_model=dataset_model,
    )

    search = af.Nautilus(
        name="mass_total[1]",
        **settings_search.search_dict,
        n_live=200 + n_live_multipole,
        n_batch=n_batch,
    )

    result = search.fit(model=model, analysis=analysis, **settings_search.fit_dict)

    return result


def run_group(
    settings_search: af.SettingsSearch,
    analysis: Union[al.AnalysisImaging, al.AnalysisInterferometer],
    source_result_for_lens: af.Result,
    source_result_for_source: af.Result,
    light_result: Optional[af.Result],
    mass_list: Optional[List] = None,
    multipole_1_list: Optional[List] = None,
    multipole_3_list: Optional[List] = None,
    multipole_4_list: Optional[List] = None,
    extra_galaxies: Optional[af.Collection] = None,
    scaling_galaxies: Optional[af.Collection] = None,
    dataset_model: Optional[af.Model] = None,
    n_batch: int = 20,
) -> af.Result:
    """
    The SLaM MASS TOTAL PIPELINE for group-scale lenses with multiple main lens galaxies.

    Group-scale counterpart to ``run``.  Fits a total mass distribution for every main
    lens simultaneously, chaining mass priors from ``source_result_for_lens`` and
    fixing light (from ``light_result``) and source (from ``source_result_for_source``).

    Parameters
    ----------
    settings_search
        AutoFit settings controlling output paths and search configuration.
    analysis
        The analysis class containing the log-likelihood function.
    source_result_for_lens
        Result whose ``model.galaxies.lens_i.mass`` and ``.shear`` provide the prior
        initialisation for each main lens mass.  Typically ``source_pix_result_1``.
    source_result_for_source
        Result whose source galaxy is fixed.  Typically ``source_pix_result_2``.
    light_result
        Result of ``light_lp.run_group``.  Provides fixed bulge/disk/point for each
        main lens.  Pass ``None`` to omit light profiles from the model.
    mass_list
        One mass model per main lens (e.g. ``af.Model(al.mp.PowerLaw)``).  Priors are
        chained from the corresponding lens in ``source_result_for_lens``.  If
        ``None``, defaults to ``af.Model(al.mp.PowerLaw)`` for every lens.
    multipole_1_list
        One optional first-order multipole per main lens.  ``None`` entries or a
        ``None`` list disable the multipole for that lens.
    multipole_3_list
        One optional third-order multipole per main lens.
    multipole_4_list
        One optional fourth-order multipole per main lens.
    extra_galaxies
        Extra galaxies with light and mass fixed from the preceding result.
    scaling_galaxies
        Scaling galaxies with light and mass fixed from the preceding result.
    dataset_model
        Optional dataset-level model components (e.g. background sky, astrometric
        offsets for multi-band fitting).
    n_batch
        Nautilus batch size.
    """
    n_lenses = sum(
        1
        for k in vars(source_result_for_lens.instance.galaxies)
        if k.startswith("lens_")
    )

    if mass_list is None:
        mass_list = [af.Model(al.mp.PowerLaw) for _ in range(n_lenses)]
    if multipole_1_list is None:
        multipole_1_list = [None] * n_lenses
    if multipole_3_list is None:
        multipole_3_list = [None] * n_lenses
    if multipole_4_list is None:
        multipole_4_list = [None] * n_lenses

    source = al.util.chaining.source_from(result=source_result_for_source)

    n_live_multipole = 0
    lens_dict = {}

    for i in range(n_lenses):
        lens_model = getattr(source_result_for_lens.model.galaxies, f"lens_{i}")

        mass = al.util.chaining.mass_from(
            mass=mass_list[i],
            mass_result=lens_model.mass,
            unfix_mass_centre=True,
        )

        multipole_1 = multipole_1_list[i]
        multipole_3 = multipole_3_list[i]
        multipole_4 = multipole_4_list[i]

        if multipole_1 is not None:
            multipole_1.m = 1
            multipole_1.centre = mass.centre
            multipole_1.einstein_radius = mass.einstein_radius
            multipole_1.slope = mass.slope
            n_live_multipole += 100

        if multipole_3 is not None:
            multipole_3.m = 3
            multipole_3.centre = mass.centre
            multipole_3.einstein_radius = mass.einstein_radius
            multipole_3.slope = mass.slope
            n_live_multipole += 100

        if multipole_4 is not None:
            multipole_4.m = 4
            multipole_4.centre = mass.centre
            multipole_4.einstein_radius = mass.einstein_radius
            multipole_4.slope = mass.slope
            n_live_multipole += 100

        if light_result is None:
            bulge = disk = point = None
        else:
            light_lens_instance = getattr(light_result.instance.galaxies, f"lens_{i}")
            bulge = light_lens_instance.bulge
            disk = light_lens_instance.disk
            point = light_lens_instance.point

        lens_dict[f"lens_{i}"] = af.Model(
            al.Galaxy,
            redshift=lens_model.redshift,
            bulge=bulge,
            disk=disk,
            point=point,
            mass=mass,
            multipole_1=multipole_1,
            multipole_3=multipole_3,
            multipole_4=multipole_4,
            shear=lens_model.shear,
        )

    lens_dict["source"] = source

    model = af.Collection(
        galaxies=af.Collection(**lens_dict),
        extra_galaxies=extra_galaxies,
        scaling_galaxies=scaling_galaxies,
        dataset_model=dataset_model,
    )

    n_live = 200 + 100 * (n_lenses - 1) + n_live_multipole

    search = af.Nautilus(
        name="mass_total[1]",
        **settings_search.search_dict,
        n_live=n_live,
        n_batch=n_batch,
    )

    return search.fit(model=model, analysis=analysis, **settings_search.fit_dict)
