import autofit as af
import autolens as al

from . import subhalo_util

from typing import Optional, Union, Tuple


def run_1_grid_search(
        settings_search: af.SettingsSearch,
        analysis: Union[al.AnalysisImaging, al.AnalysisInterferometer],
        mass_result: af.Result,
        subhalo_mass: af.Model = af.Model(al.mp.NFWMCRLudlowSph),
        free_redshift: bool = False,
        grid_dimension_arcsec: float = 3.0,
        number_of_steps: Union[Tuple[int], int] = 5,
        extra_galaxies: Optional[af.Collection] = None,
        dataset_model: Optional[af.Model] = None,
        n_batch: int = 20,
) -> af.GridSearchResult:
    """
    The SLaM SUBHALO PIPELINE for fitting lens mass models which include a dark matter subhalo.

    Parameters
    ----------
    settings_search
        The settings used to set up the non-linear search which are general to all SLaM pipelines, for example
        the `path_prefix`.
    analysis
        The analysis class which includes the `log_likelihood_function` and can be customized for the SLaM model-fit.
    mass_result
        The result of the SLaM MASS PIPELINE which ran before this pipeline.
    subhalo_result_1
        The result of the first SLaM SUBHALO PIPELINE which ran before this pipeline.
    subhalo_mass
        The `MassProfile` used to fit the subhalo in this pipeline.
    grid_dimension_arcsec
        the arc-second dimensions of the grid in the y and x directions. An input value of 3.0" means the grid in
        all four directions extends to 3.0" giving it dimensions 6.0" x 6.0".
    free_redshift
        If `True` the redshift of the subhalo is a free parameter in the second and third searches.
    number_of_steps
        The 2D dimensions of the grid (e.g. number_of_steps x number_of_steps) that the subhalo search is performed for.
    number_of_cores
        The number of cores used to perform the non-linear search grid search. If 1, each model-fit on the grid is
        performed in serial, if > 1 fits are distributed in parallel using the Python multiprocessing module.
    dataset_model
        Add aspects of the dataset to the model, for example the arc-second (y,x) offset between two datasets for
        multi-band fitting or the background sky level.
    """

    """
    __Model + Search + Analysis + Model-Fit (Search 1)__

    Search 1 of the SUBHALO PIPELINE we perform a [number_of_steps x number_of_steps] grid search of non-linear
    searches where:

     - The lens galaxy mass is modeled using MASS PIPELINE's mass distribution [Priors initialized from MASS PIPELINE].

     - The source galaxy's light is parametric or a pixelization depending on the previous MASS PIPELINE [Model and 
     priors initialized from MASS PIPELINE].

     - The subhalo redshift is fixed to that of the lens galaxy.

     - Each grid search varies the subhalo (y,x) coordinates and mass as free parameters.

     - The priors on these (y,x) coordinates are UniformPriors, with limits corresponding to the grid-cells.

    This search aims to detect a dark matter subhalo.
    """
    subhalo = af.Model(al.Galaxy, mass=subhalo_mass)

    subhalo.mass.mass_at_200 = af.LogUniformPrior(lower_limit=1.0e6, upper_limit=1.0e11)
    subhalo.mass.centre_0 = af.UniformPrior(
        lower_limit=-grid_dimension_arcsec, upper_limit=grid_dimension_arcsec
    )
    subhalo.mass.centre_1 = af.UniformPrior(
        lower_limit=-grid_dimension_arcsec, upper_limit=grid_dimension_arcsec
    )

    lens_redshift = mass_result.instance.galaxies.lens.redshift
    source_redshift = mass_result.instance.galaxies.source.redshift

    if not free_redshift:
        subhalo.redshift = lens_redshift
        subhalo.mass.redshift_object = lens_redshift
        search_tag = "search_lens_plane"
    else:
        subhalo.redshift = af.UniformPrior(lower_limit=0.0, upper_limit=source_redshift)
        subhalo.mass.redshift_object = subhalo.redshift
        search_tag = "search_multi_plane"

    subhalo.mass.redshift_source = source_redshift

    lens = mass_result.model.galaxies.lens

    source = al.util.chaining.source_from(result=mass_result)

    model = af.Collection(
        galaxies=af.Collection(lens=lens, subhalo=subhalo, source=source),
        extra_galaxies=extra_galaxies,
        dataset_model=dataset_model,
    )

    search = af.Nautilus(
        name=f"subhalo[1]_[{search_tag}]",
        **settings_search.search_dict,
        n_live=200,
        n_batch=n_batch,
    )

    subhalo_grid_search = af.SearchGridSearch(
        search=search,
        number_of_steps=number_of_steps,
    )

    result = subhalo_grid_search.fit(
        model=model,
        analysis=analysis,
        grid_priors=[
            model.galaxies.subhalo.mass.centre_1,
            model.galaxies.subhalo.mass.centre_0,
        ],
        info=settings_search.info,
    )

    try:
        subhalo_util.visualize_subhalo_detect(
            result_no_subhalo=mass_result,
            result=result,
            analysis=analysis,
            paths=subhalo_grid_search.paths,
        )
    except AttributeError:
        pass

    return result


def run_2_subhalo(
        settings_search: af.SettingsSearch,
        analysis: Union[al.AnalysisImaging, al.AnalysisInterferometer],
        mass_result: af.Result,
        subhalo_grid_search_result_1: af.GridSearchResult,
        subhalo_mass: af.Model = af.Model(al.mp.NFWMCRLudlowSph),
        free_redshift: bool = False,
        extra_galaxies: Optional[af.Collection] = None,
        dataset_model: Optional[af.Model] = None,
        n_batch: int = 20,
) -> af.Result:
    """
    The SLaM SUBHALO PIPELINE for fitting lens mass models which include a dark matter subhalo.

    Parameters
    ----------
    settings_search
        The settings used to set up the non-linear search which are general to all SLaM pipelines, for example
        the `path_prefix`.
    analysis
        The analysis class which includes the `log_likelihood_function` and can be customized for the SLaM model-fit.
    mass_result
        The result of the SLaM MASS PIPELINE which ran before this pipeline.
    subhalo_result_1
        The result of the first SLaM SUBHALO PIPELINE which ran before this pipeline.
    subhalo_mass
        The `MassProfile` used to fit the subhalo in this pipeline.
    subhalo_grid_search_result_2
        The result of the second SLaM SUBHALO PIPELINE grid search which ran before this pipeline.
    grid_dimension_arcsec
        the arc-second dimensions of the grid in the y and x directions. An input value of 3.0" means the grid in
        all four directions extends to 3.0" giving it dimensions 6.0" x 6.0".
    free_redshift
        If `True` the redshift of the subhalo is a free parameter in the second and third searches.
    dataset_model
        Add aspects of the dataset to the model, for example the arc-second (y,x) offset between two datasets for
        multi-band fitting or the background sky level.
    """

    """
    __Model + Search + Analysis + Model-Fit (Search 2)__

    Search 2 of the SUBHALO PIPELINE we refit the lens and source models above but now including a subhalo, where 
    the subhalo model is initialized from the highest evidence model of the subhalo grid search.

     - The lens galaxy mass is modeled using MASS PIPELINE's mass distribution [Priors initialized from MASS PIPELINE].

     - The source galaxy's light is parametric or a pixelization depending on the previous MASS PIPELINE [Model and 
     priors initialized from MASS PIPELINE].

     - The subhalo redshift is fixed to that of the lens galaxy.

     - Each grid search varies the subhalo (y,x) coordinates and mass as free parameters.

     - The priors on these (y,x) coordinates are UniformPriors, with limits corresponding to the grid-cells.

    This search aims to refine the parameter estimates and errors of a dark matter subhalo detected in the grid search
    above.
    """

    lens_redshift = mass_result.instance.galaxies.lens.redshift
    source_redshift = mass_result.instance.galaxies.source.redshift

    subhalo = af.Model(al.Galaxy, redshift=lens_redshift, mass=subhalo_mass)

    if not free_redshift:
        subhalo.redshift = lens_redshift
        subhalo.mass.redshift_object = lens_redshift
        refine_tag = "single_plane_refine"
    else:
        subhalo.redshift = af.UniformPrior(lower_limit=0.0, upper_limit=source_redshift)
        subhalo.mass.redshift_object = subhalo.redshift
        refine_tag = "multi_plane_refine"

    subhalo.mass.redshift_source = source_redshift
    subhalo.mass.mass_at_200 = af.LogUniformPrior(lower_limit=1.0e6, upper_limit=1.0e11)

    subhalo.mass.centre = subhalo_grid_search_result_1.model_centred_absolute(a=1.0).galaxies.subhalo.mass.centre

    # Carry through the grid-search redshift choice (important when free_redshift=True)
    subhalo.redshift = subhalo_grid_search_result_1.model.galaxies.subhalo.redshift
    subhalo.mass.redshift_object = subhalo.redshift

    model = af.Collection(
        galaxies=af.Collection(
            lens=subhalo_grid_search_result_1.model.galaxies.lens,
            subhalo=subhalo,
            source=subhalo_grid_search_result_1.model.galaxies.source,
        ),
        extra_galaxies=extra_galaxies,
        dataset_model=dataset_model,
    )

    search = af.Nautilus(
        name=f"subhalo[2]_[{refine_tag}]",
        **settings_search.search_dict,
        n_live=600,
        n_batch=n_batch,
    )

    result = search.fit(model=model, analysis=analysis, **settings_search.fit_dict)

    return result


def run_1_grid_search__multi(
        settings_search: af.SettingsSearch,
        analysis_list: Union[al.AnalysisImaging, al.AnalysisInterferometer],
        mass_result: af.Result,
        subhalo_mass: af.Model = af.Model(al.mp.NFWMCRLudlowSph),
        free_redshift: bool = False,
        grid_dimension_arcsec: float = 3.0,
        number_of_steps: Union[Tuple[int], int] = 5,
        extra_galaxies: Optional[af.Collection] = None,
        dataset_model: Optional[af.Model] = None,
        n_batch: int = 20,
) -> af.GridSearchResult:
    """
    The SLaM SUBHALO PIPELINE for fitting lens mass models which include a dark matter subhalo.

    This is a variant of the above SLaM pipeline which fits multiple datasets simultaneously using a list of
    analysis objects and the factor graph functionality. Its purpose and model are identical to the single
    dataset version, except for changes required to fit multiple datasets.

    Parameters
    ----------
    settings_search
        The settings used to set up the non-linear search which are general to all SLaM pipelines, for example
        the `path_prefix`.
    analysis
        The analysis class which includes the `log_likelihood_function` and can be customized for the SLaM model-fit.
    mass_result
        The result of the SLaM MASS PIPELINE which ran before this pipeline.
    subhalo_result_1
        The result of the first SLaM SUBHALO PIPELINE which ran before this pipeline.
    subhalo_mass
        The `MassProfile` used to fit the subhalo in this pipeline.
    grid_dimension_arcsec
        the arc-second dimensions of the grid in the y and x directions. An input value of 3.0" means the grid in
        all four directions extends to 3.0" giving it dimensions 6.0" x 6.0".
    free_redshift
        If `True` the redshift of the subhalo is a free parameter in the second and third searches.
    number_of_steps
        The 2D dimensions of the grid (e.g. number_of_steps x number_of_steps) that the subhalo search is performed for.
    number_of_cores
        The number of cores used to perform the non-linear search grid search. If 1, each model-fit on the grid is
        performed in serial, if > 1 fits are distributed in parallel using the Python multiprocessing module.
    dataset_model
        Add aspects of the dataset to the model, for example the arc-second (y,x) offset between two datasets for
        multi-band fitting or the background sky level.
    """

    """
    __Model + Search + Analysis + Model-Fit (Search 1)__

    Search 1 of the SUBHALO PIPELINE we perform a [number_of_steps x number_of_steps] grid search of non-linear
    searches where:

     - The lens galaxy mass is modeled using MASS PIPELINE's mass distribution [Priors initialized from MASS PIPELINE].

     - The source galaxy's light is parametric or a pixelization depending on the previous MASS PIPELINE [Model and 
     priors initialized from MASS PIPELINE].

     - The subhalo redshift is fixed to that of the lens galaxy.

     - Each grid search varies the subhalo (y,x) coordinates and mass as free parameters.

     - The priors on these (y,x) coordinates are UniformPriors, with limits corresponding to the grid-cells.

    This search aims to detect a dark matter subhalo.
    """

    subhalo = af.Model(al.Galaxy, mass=subhalo_mass)

    subhalo.mass.mass_at_200 = af.LogUniformPrior(lower_limit=1.0e6, upper_limit=1.0e11)
    subhalo.mass.centre_0 = af.UniformPrior(
        lower_limit=-grid_dimension_arcsec, upper_limit=grid_dimension_arcsec
    )
    subhalo.mass.centre_1 = af.UniformPrior(
        lower_limit=-grid_dimension_arcsec, upper_limit=grid_dimension_arcsec
    )

    lens_redshift = mass_result[0].instance.galaxies.lens.redshift
    source_redshift = mass_result[0].instance.galaxies.source.redshift

    if not free_redshift:
        subhalo.redshift = lens_redshift
        subhalo.mass.redshift_object = lens_redshift
        search_tag = "search_lens_plane"
    else:
        subhalo.redshift = af.UniformPrior(lower_limit=0.0, upper_limit=source_redshift)
        subhalo.mass.redshift_object = subhalo.redshift
        search_tag = "search_multi_plane"

    subhalo.mass.redshift_source = source_redshift

    analysis_factor_list = []

    for i, analysis in enumerate(analysis_list):
        lens = mass_result[i].model.galaxies.lens
        source = al.util.chaining.source_from(result=mass_result[i])

        model = af.Collection(
            galaxies=af.Collection(lens=lens, subhalo=subhalo, source=source),
            extra_galaxies=extra_galaxies,
            dataset_model=dataset_model,
        )

        analysis_factor_list.append(af.AnalysisFactor(prior_model=model, analysis=analysis))

    factor_graph = af.FactorGraphModel(*analysis_factor_list, use_jax=True)

    search = af.Nautilus(
        name=f"subhalo[1]_[{search_tag}]",
        **settings_search.search_dict,
        n_live=200,
        n_batch=n_batch,
    )

    subhalo_grid_search = af.SearchGridSearch(
        search=search,
        number_of_steps=number_of_steps,
    )

    # GridSearch currently expects a (single) prior model to define the grid priors;
    # use the factor graph's global prior model for that.
    global_model = factor_graph.global_prior_model

    result = subhalo_grid_search.fit(
        model=global_model,
        analysis=factor_graph,
        grid_priors=[
            global_model.galaxies.subhalo.mass.centre_1,
            global_model.galaxies.subhalo.mass.centre_0,
        ],
        info=settings_search.info,
    )

    try:
        subhalo_util.visualize_subhalo_detect(
            result_no_subhalo=mass_result[0],
            result=result,
            analysis=analysis_list[0],
            paths=subhalo_grid_search.paths,
        )
    except AttributeError:
        pass

    return result


def run_2_subhalo__multi(
        settings_search: af.SettingsSearch,
        analysis_list: Union[al.AnalysisImaging, al.AnalysisInterferometer],
        mass_result: af.Result,
        subhalo_grid_search_result_1: af.GridSearchResult,
        subhalo_mass: af.Model = af.Model(al.mp.NFWMCRLudlowSph),
        free_redshift: bool = False,
        extra_galaxies: Optional[af.Collection] = None,
        dataset_model: Optional[af.Model] = None,
        n_batch: int = 20,
) -> af.Result:
    """
    The SLaM SUBHALO PIPELINE for fitting lens mass models which include a dark matter subhalo.

    This is a variant of the above SLaM pipeline which fits multiple datasets simultaneously using a list of
    analysis objects and the factor graph functionality. Its purpose and model are identical to the single
    dataset version, except for changes required to fit multiple datasets.

    Parameters
    ----------
    settings_search
        The settings used to set up the non-linear search which are general to all SLaM pipelines, for example
        the `path_prefix`.
    analysis
        The analysis class which includes the `log_likelihood_function` and can be customized for the SLaM model-fit.
    mass_result
        The result of the SLaM MASS PIPELINE which ran before this pipeline.
    subhalo_result_1
        The result of the first SLaM SUBHALO PIPELINE which ran before this pipeline.
    subhalo_mass
        The `MassProfile` used to fit the subhalo in this pipeline.
    subhalo_grid_search_result_2
        The result of the second SLaM SUBHALO PIPELINE grid search which ran before this pipeline.
    grid_dimension_arcsec
        the arc-second dimensions of the grid in the y and x directions. An input value of 3.0" means the grid in
        all four directions extends to 3.0" giving it dimensions 6.0" x 6.0".
    free_redshift
        If `True` the redshift of the subhalo is a free parameter in the second and third searches.
    dataset_model
        Add aspects of the dataset to the model, for example the arc-second (y,x) offset between two datasets for
        multi-band fitting or the background sky level.
    """

    """
    __Model + Search + Analysis + Model-Fit (Search 2)__

    Search 2 of the SUBHALO PIPELINE we refit the lens and source models above but now including a subhalo, where 
    the subhalo model is initialized from the highest evidence model of the subhalo grid search.

     - The lens galaxy mass is modeled using MASS PIPELINE's mass distribution [Priors initialized from MASS PIPELINE].

     - The source galaxy's light is parametric or a pixelization depending on the previous MASS PIPELINE [Model and 
     priors initialized from MASS PIPELINE].

     - The subhalo redshift is fixed to that of the lens galaxy.

     - Each grid search varies the subhalo (y,x) coordinates and mass as free parameters.

     - The priors on these (y,x) coordinates are UniformPriors, with limits corresponding to the grid-cells.

    This search aims to refine the parameter estimates and errors of a dark matter subhalo detected in the grid search
    above.
    """

    lens_redshift = mass_result[0].instance.galaxies.lens.redshift
    source_redshift = mass_result[0].instance.galaxies.source.redshift

    subhalo = af.Model(al.Galaxy, redshift=lens_redshift, mass=subhalo_mass)

    if not free_redshift:
        subhalo.redshift = lens_redshift
        subhalo.mass.redshift_object = lens_redshift
        refine_tag = "single_plane_refine"
    else:
        subhalo.redshift = af.UniformPrior(lower_limit=0.0, upper_limit=source_redshift)
        subhalo.mass.redshift_object = subhalo.redshift
        refine_tag = "multi_plane_refine"

    subhalo.mass.redshift_source = source_redshift
    subhalo.mass.mass_at_200 = af.LogUniformPrior(lower_limit=1.0e6, upper_limit=1.0e11)

    subhalo.mass.centre = subhalo_grid_search_result_1.model_centred_absolute(a=1.0).galaxies.subhalo.mass.centre

    subhalo.redshift = subhalo_grid_search_result_1.model.galaxies.subhalo.redshift
    subhalo.mass.redshift_object = subhalo.redshift

    analysis_factor_list = []

    for analysis in analysis_list:
        model = af.Collection(
            galaxies=af.Collection(
                lens=subhalo_grid_search_result_1.model.galaxies.lens,
                subhalo=subhalo,
                source=subhalo_grid_search_result_1.model.galaxies.source,
            ),
            extra_galaxies=extra_galaxies,
            dataset_model=dataset_model,
        )

        analysis_factor_list.append(af.AnalysisFactor(prior_model=model, analysis=analysis))

    factor_graph = af.FactorGraphModel(*analysis_factor_list, use_jax=True)

    search = af.Nautilus(
        name=f"subhalo[2]_[{refine_tag}]",
        **settings_search.search_dict,
        n_live=600,
        n_batch=n_batch,
    )

    result = search.fit(model=factor_graph.global_prior_model, analysis=factor_graph)

    return result
