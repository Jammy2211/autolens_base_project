import autofit as af
import autolens as al

from typing import Optional, Tuple, Union


def run_1(
    settings_search: af.SettingsSearch,
    analysis: Union[al.AnalysisImaging, al.AnalysisInterferometer],
    source_lp_result: af.Result,
    mesh_init: af.Model(al.AbstractMesh) = af.Model(al.mesh.Delaunay),
    regularization_init: af.Model(al.AbstractRegularization) = af.Model(
        al.reg.AdaptSplit
    ),
    extra_galaxies: Optional[af.Collection] = None,
    dataset_model: Optional[af.Model] = None,
    fixed_mass_model: bool = False,
    n_batch: int = 20,
) -> af.Result:
    """
    The first SLaM SOURCE PIX PIPELINE, which initializes a lens model which uses a pixelized source for the source
    analysis.

    The first SOURCE PIX PIPELINE may require an adapt-image, for example to adapt the regularization scheme to the
    source's unlensed morphology. The adapt image provided by the SOURCE LP PIPELINE may not cover the entire source
    galaxy (e.g. because the MGE only captures part of the source) and produce a suboptimal fit.

    The result of this pipeline is used in the second SOURCE PIX PIPELINE to adapt the source pixelization to the
    source's unlensed morphology via an adapt image, where the adapt image produced in this pipeline will give a robust
    source image because it uses a pixelized source.

    Parameters
    ----------
    settings_search
        The settings used to set up the non-linear search which are general to all SLaM pipelines, for example
        the `path_prefix`.
    analysis
        The analysis class which includes the `log_likelihood_function` and can be customized for the SLaM model-fit.
    source_lp_result
        The results of the SLaM SOURCE LP PIPELINE which ran before this pipeline.
    mesh_init
        The mesh, which defines how the source is reconstruction in the source-plane, used by the pixelization
        in the first search which initializes the source.
    regularization_init
        The regularization, which places a smoothness prior on the source reconstruction, used by the pixelization
        which fits the source light in the initialization search (`search[1]`).
    extra_galaxies
        Additional extra galaxies containing light and mass profiles, which model nearby line of sight galaxies.
    dataset_model
        Add aspects of the dataset to the model, for example the arc-second (y,x) offset between two datasets for
        multi-band fitting or the background sky level.
    fixed_mass_model
        Whether the mass model is fixed from the SOURCE LP PIPELINE, which is generally used for multi-band fitting
        where the mass model is fixed to the first band, albeit it may work for standard fitting if the SOURCE LP
        PIPELINE provides a good mass model.
    """

    """
    __Model + Search + Analysis + Model-Fit (Search 1)__

    Search 1 of the SOURCE PIX PIPELINE fits a lens model where:

    - The lens galaxy light is modeled using light profiles [parameters fixed to result of SOURCE LP PIPELINE].

     - The lens galaxy mass is modeled using a total mass distribution [model initialized from the results of the 
     SOURCE LP PIPELINE].

     - The source galaxy's light is the input initialization image mesh, mesh and regularization scheme [parameters of 
     regularization free to vary].

    This search improves the lens mass model by modeling the source using a pixelization and computes the adapt
    images that are used in search 2.
    """

    if not fixed_mass_model:
        mass = al.util.chaining.mass_from(
            mass=source_lp_result.model.galaxies.lens.mass,
            mass_result=source_lp_result.model.galaxies.lens.mass,
            unfix_mass_centre=True,
        )
        shear = source_lp_result.model.galaxies.lens.shear

    else:
        mass = source_lp_result.instance.galaxies.lens.mass
        shear = source_lp_result.instance.galaxies.lens.shear

    model = af.Collection(
        galaxies=af.Collection(
            lens=af.Model(
                al.Galaxy,
                redshift=source_lp_result.instance.galaxies.lens.redshift,
                bulge=source_lp_result.instance.galaxies.lens.bulge,
                disk=source_lp_result.instance.galaxies.lens.disk,
                point=source_lp_result.instance.galaxies.lens.point,
                mass=mass,
                shear=shear,
            ),
            source=af.Model(
                al.Galaxy,
                redshift=source_lp_result.instance.galaxies.source.redshift,
                pixelization=af.Model(
                    al.Pixelization,
                    mesh=mesh_init,
                    regularization=regularization_init,
                ),
            ),
        ),
        extra_galaxies=extra_galaxies,
        dataset_model=dataset_model,
    )

    search = af.Nautilus(
        name="source_pix[1]",
        **settings_search.search_dict,
        n_live=150,
        n_batch=n_batch,
    )

    result = search.fit(model=model, analysis=analysis, **settings_search.fit_dict)

    return result


def run_1_group(
    settings_search: af.SettingsSearch,
    analysis: Union[al.AnalysisImaging, al.AnalysisInterferometer],
    source_lp_result: af.Result,
    mesh_init: af.Model(al.AbstractMesh) = af.Model(al.mesh.Delaunay),
    regularization_init: af.Model(al.AbstractRegularization) = af.Model(
        al.reg.AdaptSplit
    ),
    extra_galaxies: Optional[af.Collection] = None,
    scaling_galaxies: Optional[af.Collection] = None,
    dataset_model: Optional[af.Model] = None,
    fixed_mass_model: bool = False,
    n_batch: int = 20,
) -> af.Result:
    """
    The first SLaM SOURCE PIX PIPELINE for group-scale lenses with multiple main lens
    galaxies.

    Group-scale counterpart to ``run_1``.  The number of main lenses is inferred from
    the result so the function works for any number of lenses.

    For each main lens:
      - bulge / disk / point fixed to SOURCE LP result instances.
      - mass and shear free, chained from SOURCE LP result model
        (or fixed if ``fixed_mass_model=True``).

    Parameters
    ----------
    settings_search
        AutoFit settings controlling output paths and search configuration.
    analysis
        The analysis class containing the log-likelihood function.
    source_lp_result
        Result of ``source_lp.run_group``.  Must contain ``galaxies.lens_0``, etc.
    mesh_init
        Mesh used by the pixelization in this initialization search.
    regularization_init
        Regularization scheme used in this initialization search.
    extra_galaxies
        Combined extra + scaling galaxy collection passed through from the caller.
    dataset_model
        Optional dataset-level model components.
    fixed_mass_model
        If ``True``, fix mass and shear to the SOURCE LP instance rather than
        treating them as free models.
    """
    n_lenses = sum(
        1 for k in vars(source_lp_result.instance.galaxies) if k.startswith("lens_")
    )

    lens_dict = {}
    for i in range(n_lenses):
        lp_lens_instance = getattr(source_lp_result.instance.galaxies, f"lens_{i}")
        lp_lens_model = getattr(source_lp_result.model.galaxies, f"lens_{i}")

        if not fixed_mass_model:
            mass = al.util.chaining.mass_from(
                mass=lp_lens_model.mass,
                mass_result=lp_lens_model.mass,
                unfix_mass_centre=True,
            )
            shear = lp_lens_model.shear
        else:
            mass = lp_lens_instance.mass
            shear = lp_lens_instance.shear

        lens_dict[f"lens_{i}"] = af.Model(
            al.Galaxy,
            redshift=lp_lens_instance.redshift,
            bulge=lp_lens_instance.bulge,
            disk=lp_lens_instance.disk,
            point=lp_lens_instance.point,
            mass=mass,
            shear=shear,
        )

    lens_dict["source"] = af.Model(
        al.Galaxy,
        redshift=source_lp_result.instance.galaxies.source.redshift,
        pixelization=af.Model(
            al.Pixelization,
            mesh=mesh_init,
            regularization=regularization_init,
        ),
    )

    model = af.Collection(
        galaxies=af.Collection(**lens_dict),
        extra_galaxies=extra_galaxies,
        scaling_galaxies=scaling_galaxies,
        dataset_model=dataset_model,
    )

    # n_live scales with the number of free mass models.
    n_live = 150 + 50 * (n_lenses - 1)

    search = af.Nautilus(
        name="source_pix[1]",
        **settings_search.search_dict,
        n_live=n_live,
        n_batch=n_batch,
    )

    result = search.fit(model=model, analysis=analysis, **settings_search.fit_dict)

    return result


def run_2_group(
    settings_search: af.SettingsSearch,
    analysis: Union[al.AnalysisImaging, al.AnalysisInterferometer],
    source_lp_result: af.Result,
    source_pix_result_1: af.Result,
    mesh: af.Model(al.AbstractMesh) = af.Model(al.mesh.Delaunay),
    regularization: af.Model(al.AbstractRegularization) = af.Model(al.reg.AdaptSplit),
    dataset_model: Optional[af.Model] = None,
    n_batch: int = 20,
) -> af.Result:
    """
    The second SLaM SOURCE PIX PIPELINE for group-scale lenses with multiple main lens
    galaxies.

    Group-scale counterpart to ``run_2``.  All lens components are fixed; only the
    pixelization mesh and regularization are free.

    For each main lens:
      - bulge / disk / point fixed to SOURCE LP result instances.
      - mass and shear fixed to SOURCE PIX[1] result instances.

    ``extra_galaxies`` are taken directly from ``source_pix_result_1.instance``,
    matching the behaviour of the single-lens ``run_2``.

    Parameters
    ----------
    settings_search
        AutoFit settings controlling output paths and search configuration.
    analysis
        The analysis class containing the log-likelihood function.
    source_lp_result
        Result of ``source_lp.run_group``.  Provides light and source redshift.
    source_pix_result_1
        Result of ``run_1_group``.  Provides fixed mass/shear and extra_galaxies.
    mesh
        Final mesh for the pixelization.
    regularization
        Final regularization scheme.
    dataset_model
        Optional dataset-level model components.
    """
    n_lenses = sum(
        1 for k in vars(source_lp_result.instance.galaxies) if k.startswith("lens_")
    )

    lens_dict = {}
    for i in range(n_lenses):
        lp_lens_instance = getattr(source_lp_result.instance.galaxies, f"lens_{i}")
        pix1_lens_instance = getattr(source_pix_result_1.instance.galaxies, f"lens_{i}")

        lens_dict[f"lens_{i}"] = af.Model(
            al.Galaxy,
            redshift=lp_lens_instance.redshift,
            bulge=lp_lens_instance.bulge,
            disk=lp_lens_instance.disk,
            point=lp_lens_instance.point,
            mass=pix1_lens_instance.mass,
            shear=pix1_lens_instance.shear,
        )

    lens_dict["source"] = af.Model(
        al.Galaxy,
        redshift=source_lp_result.instance.galaxies.source.redshift,
        pixelization=af.Model(
            al.Pixelization,
            mesh=mesh,
            regularization=regularization,
        ),
    )

    model = af.Collection(
        galaxies=af.Collection(**lens_dict),
        extra_galaxies=source_pix_result_1.instance.extra_galaxies,
        scaling_galaxies=source_pix_result_1.instance.scaling_galaxies,
        dataset_model=dataset_model,
    )

    search = af.Nautilus(
        name="source_pix[2]",
        **settings_search.search_dict,
        n_live=75,
        n_batch=n_batch,
    )

    result = search.fit(model=model, analysis=analysis, **settings_search.fit_dict)

    return result


def run_2(
    settings_search: af.SettingsSearch,
    analysis: Union[al.AnalysisImaging, al.AnalysisInterferometer],
    source_lp_result: af.Result,
    source_pix_result_1: af.Result,
    mesh: af.Model(al.AbstractMesh) = af.Model(al.mesh.Delaunay),
    regularization: af.Model(al.AbstractRegularization) = af.Model(al.reg.AdaptSplit),
    dataset_model: Optional[af.Model] = None,
    n_batch: int = 20,
) -> af.Result:
    """
    The second SLaM SOURCE PIX PIPELINE, which fits a fixed lens model which uses a pixelized source for the source
    analysis.

    The second SOURCE PIX PIPELINE performs a fit using an advanced pixelizaiton which adapt the source's pixelization
    to the source's unlensed morphology.

    This feature requires an adapt-image, which is computed after the first SOURCE PIX PIPELINE.

    Parameters
    ----------
    settings_search
        The settings used to set up the non-linear search which are general to all SLaM pipelines, for example
        the `path_prefix`.
    analysis
        The analysis class which includes the `log_likelihood_function` and can be customized for the SLaM model-fit.
    source_lp_result
        The results of the SLaM SOURCE LP PIPELINE which ran before this pipeline.
    image_mesh
        The image mesh, which defines how the mesh centres are computed in the image-plane, used by the pixelization
        in the final search which improves the source adaption.
    mesh
        The mesh, which defines how the source is reconstruction in the source-plane, used by the pixelization
        in the final search which improves the source adaption.
    regularization
        The regularization, which places a smoothness prior on the source reconstruction, used by the pixelization
        in the final search which improves the source adaption.
    extra_galaxies
        Additional extra galaxies containing light and mass profiles, which model nearby line of sight galaxies.
    dataset_model
        Add aspects of the dataset to the model, for example the arc-second (y,x) offset between two datasets for
        multi-band fitting or the background sky level.
    """

    """
    __Model + Search + Analysis + Model-Fit (Search 2)__

    Search 2 of the SOURCE PIX PIPELINE fits a lens model where:

    - The lens galaxy light is modeled using a light profiles [parameters fixed to result of SOURCE LP PIPELINE].
    - The lens galaxy mass is modeled using a total mass distribution [parameters fixed to result of search 1].
    - The source galaxy's light is the input final mesh and regularization.

    This search initializes the pixelization's mesh and regularization.
    """
    model = af.Collection(
        galaxies=af.Collection(
            lens=af.Model(
                al.Galaxy,
                redshift=source_lp_result.instance.galaxies.lens.redshift,
                bulge=source_lp_result.instance.galaxies.lens.bulge,
                disk=source_lp_result.instance.galaxies.lens.disk,
                point=source_lp_result.instance.galaxies.lens.point,
                mass=source_pix_result_1.instance.galaxies.lens.mass,
                shear=source_pix_result_1.instance.galaxies.lens.shear,
            ),
            source=af.Model(
                al.Galaxy,
                redshift=source_lp_result.instance.galaxies.source.redshift,
                pixelization=af.Model(
                    al.Pixelization,
                    mesh=mesh,
                    regularization=regularization,
                ),
            ),
        ),
        extra_galaxies=source_pix_result_1.instance.extra_galaxies,
        dataset_model=dataset_model,
    )

    search = af.Nautilus(
        name="source_pix[2]",
        **settings_search.search_dict,
        n_live=75,
        n_batch=n_batch,
    )

    result = search.fit(model=model, analysis=analysis, **settings_search.fit_dict)

    return result


def run_1__bypass_lp(
    settings_search: af.SettingsSearch,
    analysis: Union[al.AnalysisImaging, al.AnalysisInterferometer],
    lens_bulge: Optional[af.Model] = af.Model(al.lp.Sersic),
    lens_disk: Optional[af.Model] = af.Model(al.lp.Exponential),
    lens_point: Optional[af.Model] = None,
    mass: af.Model = af.Model(al.mp.Isothermal),
    shear: af.Model(al.mp.ExternalShear) = af.Model(al.mp.ExternalShear),
    mesh_init: af.Model(al.AbstractMesh) = af.Model(al.mesh.Delaunay),
    regularization_init: af.Model(al.AbstractRegularization) = af.Model(
        al.reg.AdaptSplit
    ),
    redshift_lens: float = 0.5,
    redshift_source: float = 1.0,
    extra_galaxies: Optional[af.Collection] = None,
    dataset_model: Optional[af.Model] = None,
    n_batch: int = 20,
) -> af.Result:
    """
    The first SLaM SOURCE PIX PIPELINE, which initializes a lens model which uses a pixelized source for the source
    analysis.

    This variant bypasses the source lp pipeline and is typically used for the interferometer SLaM pipeline, albeit
    it can also be used for imaging analysis.

    The first SOURCE PIX PIPELINE may require an adapt-image, for example to adapt the regularization scheme to the
    source's unlensed morphology. The adapt image provided by the SOURCE LP PIPELINE may not cover the entire source
    galaxy (e.g. because the MGE only captures part of the source) and produce a suboptimal fit.

    The result of this pipeline is used in the second SOURCE PIX PIPELINE to adapt the source pixelization to the
    source's unlensed morphology via an adapt image, where the adapt image produced in this pipeline will give a robust
    source image because it uses a pixelized source.

    Parameters
    ----------
    settings_search
        The settings used to set up the non-linear search which are general to all SLaM pipelines, for example
        the `path_prefix`.
    analysis
        The analysis class which includes the `log_likelihood_function` and can be customized for the SLaM model-fit.
    mesh_init
        The mesh, which defines how the source is reconstruction in the source-plane, used by the pixelization
        in the first search which initializes the source.
    regularization_init
        The regularization, which places a smoothness prior on the source reconstruction, used by the pixelization
        which fits the source light in the initialization search (`search[1]`).
    extra_galaxies
        Additional extra galaxies containing light and mass profiles, which model nearby line of sight galaxies.
    dataset_model
        Add aspects of the dataset to the model, for example the arc-second (y,x) offset between two datasets for
        multi-band fitting or the background sky level.
    fixed_mass_model
        Whether the mass model is fixed from the SOURCE LP PIPELINE, which is generally used for multi-band fitting
        where the mass model is fixed to the first band, albeit it may work for standard fitting if the SOURCE LP
        PIPELINE provides a good mass model.
    """

    """
    __Model + Search + Analysis + Model-Fit (Search 1)__

    Search 1 of the SOURCE PIX PIPELINE fits a lens model where:

    - The lens galaxy light is modeled using light profiles [parameters fixed to result of SOURCE LP PIPELINE].

     - The lens galaxy mass is modeled using a total mass distribution [model initialized from the results of the 
     SOURCE LP PIPELINE].

     - The source galaxy's light is the input initialization image mesh, mesh and regularization scheme [parameters of 
     regularization free to vary].

    This search improves the lens mass model by modeling the source using a pixelization and computes the adapt
    images that are used in search 2.
    """

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
                pixelization=af.Model(
                    al.Pixelization,
                    mesh=mesh_init,
                    regularization=regularization_init,
                ),
            ),
        ),
        extra_galaxies=extra_galaxies,
        dataset_model=dataset_model,
    )

    search = af.Nautilus(
        name="source_pix[1]",
        **settings_search.search_dict,
        n_live=150,
        n_batch=n_batch,
    )

    result = search.fit(model=model, analysis=analysis, **settings_search.fit_dict)

    return result

