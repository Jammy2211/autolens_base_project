---
title: "`PyAutoLens-Assistant`: Using Natural Language and AI to Analyse Gravitational Lenses"
tags:
  - Python
  - astronomy
  - gravitational lensing
  - artificial intelligence
  - agentic ai
  - large language models
  - natural language interfaces
authors:
  - name: James W. Nightingale
    orcid: 0000-0002-8987-7401
    affiliation: 1
    corresponding: true
    email: James.Nightingale@Newcastle.ac.uk
  - name: Richard G. Hayes
    affiliation: 2
affiliations:
  - name: School of Mathematics, Statistics and Physics, Newcastle University, Herschel Building, Newcastle upon Tyne, NE1 7RU, United Kingdom
    index: 1
  - name: Institute for Computational Cosmology, Durham University, South Road, Durham DH1 3LE, United Kingdom
    index: 2
date: 30 July 2026
submitted_at: 2026-07-30
bibliography: paper.bib
software_repository_url: https://github.com/PyAutoLabs/autolens_assistant
---

# Summary

Stage IV weak-lensing surveys, such as Euclid [@EuclidCollaboration2025] and the Vera C. Rubin 
Observatory [@LSSTDarkEnergyScienceCollaboration2012], are mapping the distribution of mass across the Universe on 
an unprecedented scale. They include strong-lens searches which are rapidly expanding samples of galaxy-, group-, and 
cluster-scale lenses beyond hundreds of thousands [@Collett2015]. Lensing studies draw on optical and infrared imaging [@Bolton2006; @Nightingale2025COWLS], 
submm and radio interferometry [@Hezaveh2016; @Vegetti2025; @Rizzo2020], strongly lensed point sources and 
transients such as quasars and supernovae [@Wong2019; @Grillo2018], and weak-lensing shear catalogues. Together, these datasets enable studies of cosmology, dark matter, galaxy formation, 
star formation, and the early Universe. [`PyAutoLens-JAX`](https://github.com/PyAutoLabs/PyAutoLens) [@Nightingale2021; @NightingaleJAX2026] provides open-source software for GPU-native, 
autodifferentiable joint lensing analyses across these datasets and scales. However, combining datasets and performing 
inference across vastly different lensing scales is inherently complex, requiring 
substantial effort from the scientist to find and adapt the `PyAutoLens-JAX` API and Python syntax to build their specific 
analysis pipelines.

[`PyAutoLens-Assistant`](https://github.com/PyAutoLabs/autolens_assistant) allows scientists to use natural language to learn and perform
strong- and weak-lensing analyses across galaxy, group, and cluster scales. It runs inside an AI coding agent such as 
`Claude Code` or `Codex`: users can ask `PyAutoLens-Assistant` how a lensing analysis works, plan a study or learn the 
physics in its teacher mode without executing anything, and obtain fully documented code for the lens analysis task; 
when they describe the desired analysis, the agent collates the data, writes and executes Python scripts and brings 
together the results. The scientist can then use natural language to visualize, investigate and interpret the results. `PyAutoLens-Assistant` also 
includes a Teacher Mode for users new to `PyAutoLens-JAX`  or gravitational lensing, which can explain the core domain specific 
concepts whilest directing them to documentation and Jupyter Notebook guides to help them further build their understanding.

# Statement of need

Scientists can often inspect a dataset and know exactly the analysis they want to perform. For example, 
when shown the multi-wavelength James Webb Space Telescope (JWST) observations of the 
COSMOS-Web Ring [@Casey2023; @Mercier2024] in \autoref{fig:cosmos_web_ring_dataset}:

![The COSMOS-Web Ring dataset. The left panel is an RGB colour composite from COSMOS-Web, annotated with the 
lens galaxy, the lensed source galaxy and a nearby extra galaxy. The middle and right panels show the JWST 
F277W and F444W imaging, with the dashed circle marking the $1.8''$ mask applied in the analysis. 
\label{fig:cosmos_web_ring_dataset}](cosmos_web_ring_dataset.png)

An experienced strong-lensing scientist might examine this data and say:

> I want to model the F277W and F444W JWST imaging of the COSMOS-Web Ring independently, which are in 
> my local folder dataset/cosmos_web_ring. Model the lens light with a multi-Gaussian expansion, its mass with a 
> singular isothermal ellipsoid plus external shear, and model the source using a multi-Gaussian expansion. For 
> speed, run the analysis on the laptop GPU using a JAX optimizer that estimates only the maximum-likelihood solution. 
> Plot the observed image at each wavelength in the left column, its lensed source model in the middle column, and its source 
> reconstruction in the right column.

When the prompt above is input into `PyAutoLens-Assistant` using `Claude Code` Opus 5.0, after answering a couple of 
clarifying questions, the end-to-end analysis produces \autoref{fig:cosmos_web_ring}, successfully delivering the 
output requested in the prompt.

The example above is a **natural-language workflow**: the scientist specifies the analysis in scientific terms, and the 
assistant translates it into executable code. Complex modelling concepts can therefore be expressed clearly even 
when implementing them manually in code would require substantial effort. By handling the Python syntax and `PyAutoLens-JAX` API, 
the assistant allows scientists to focus on what the analysis should do and why, rather than facing the implementation burden of 
how to implement it in code. Its GitHub page has several example prompts which illustrate this further. One uses a 
three-paragraph prompt to reproduce the well-known detection of a dark matter subhalo in the strong lens SDSS J0946+1006 [@Vegetti2010; @Minor2021]; 
another simulates CCD imaging, interferometric and weak lensing observations of a multi-galaxy strong lens and then 
performs joint inference of all datasets.

As lens analysis becomes more complex, it is equally important that performing it becomes **more accessible**.
With `PyAutoLens-Assistant`, a user does not need prior expertise in gravitational lensing.
For example, the following prompt initiates a guided conversation in which the assistant explains the available 
analysis approaches step by step:


> Find the data for the COSMOS-Web Ring, provide a short `PyAutoLens-JAX` script to plot it, and, given that I am a new 
> user, explain the different ways in which we could model this strong-lensing system.

![The end-to-end COSMOS-Web Ring analysis produced by `PyAutoLens-Assistant` from the natural-language prompt above. 
Each row corresponds to one JWST band (F277W and F444W); the left column shows the observed image, the 
middle column the lensed source model, and the right column the source-plane reconstruction. 
\label{fig:cosmos_web_ring}](model_fit.png)

The above prompt will activate `PyAutoLens-Assistant`'s Teacher Mode, which explains scientific concepts in context and directs users to relevant 
documentation, tutorials, and open-source lectures. It can introduce undergraduate and early-stage PhD researchers 
to scientific data analysis and Bayesian inference, while helping experienced lensing scientists who are new 
to `PyAutoLens-JAX` navigate its API and workflows. This support is increasingly valuable as lensing datasets span 
multiple observing techniques, physical scales, and strong- and weak-lensing regimes. By reducing the need to work 
directly with software syntax, `PyAutoLens-Assistant` enables users to learn through real analyses while focusing on 
the underlying scientific concepts.

# How it works

`PyAutoLens-Assistant` builds on the general knowledge and reasoning capabilities of its underlying foundation model by 
adding two version-controlled packages containing `PyAutoLens-JAX`-specific context. These comprise libraries of Markdown files 
structured for rapid machine reading, enabling the assistant to consult the `PyAutoLens-JAX` documentation,
end-to-end examples, and relevant scientific literature when responding to a user.

The first package is a reference wiki comprising core and literature pages. The core wiki provides an AI-readable 
interface to the `PyAutoLens-JAX` API documentation, so that the assistant can identify and retrieve the information 
and required API relevant to a user’s request. The literature wiki is derived from more than 300 strong-lensing 
papers and summarises scientific concepts, terminology, surveys, and studies of individual lens systems. This allows 
the assistant to interpret an input prompt within the wider scientific literature. Users can add papers relevant to 
their own research, giving the assistant a bespoke knowledge base tailored to their particular scientific study.

The second package is the skills library, which provides the procedural knowledge required to perform end-to-end analyses. 
A skill is a structured instruction file used by AI coding agents to complete a specific task. During more than a 
decade of `PyAutoLens-JAX` development, over 300 human-written examples have been created 
in the [`autolens_workspace`](https://github.com/PyAutoLabs/autolens_workspace), covering datasets such as CCD imaging and interferometric observations and tasks including 
simulation, lens modelling and result analysis. These examples are extensively documented, containing 
substantially more explanatory text than code, and therefore describing both the implementation and the scientific 
reasoning behind each step. `PyAutoLens-Assistant` distils this material into concise skill files that describe how to 
perform specific tasks, which direct the assistant to the full workspace examples if more context is required.

`PyAutoLens-Assistant` does not simply match a user’s prompt to a single workspace example. The example prompts on the 
`PyAutoLens-Assistant` GitHub page deliberately request analyses that no individual `autolens_workspace` example performs. Instead, the assistant 
combines the relevant skills and reference material, generalizing across multiple examples to construct an end-to-end 
workflow tailored to the user’s dataset and scientific aims. This design also makes the assistant straightforward to 
maintain: as new `PyAutoLens-JAX` features and accompanying `autolens_workspace` examples are introduced, the corresponding 
skills are updated so that the assistant immediately supports the new functionality in a way that generalizes
across all other functionality.

`PyAutoLens-Assistant` can be used through an AI chat assistant such as `ChatGPT` or `Claude`. 
The [`llms.txt`](https://github.com/PyAutoLabs/autolens_assistant/blob/main/llms.txt) file defines the canonical reading order through the wikis, skills, and runnable examples, enabling 
the assistant to answer questions, explain scientific concepts, locate relevant documentation, and generate end-to-end 
modelling scripts. This mode is necessarily manual and limited: the assistant cannot normally inspect local files or 
execute code, so users must run scripts themselves and copy code, errors, figures, and other outputs between the 
conversation and their working environment. Nevertheless, conversational AI is currently the interface most familiar 
to astronomers, making this mode easy to integrate into existing workflows.

AI CLI coding agents such as `Claude Code` and `Codex` provide the complete `PyAutoLens-Assistant` workflow using the same 
curated knowledge base with substantially fewer manual steps. They can inspect `.fits` datasets directly (e.g. reading
their header information), write and execute end-to-end analysis scripts, read and then resolve Python exceptions
if code runs incorrectly. `PyAutoLens-JAX` outputs lens modeling results into directories locally which include 
structured `.json` files describing the model, sampler, results and other metadata which are specifically designed to 
be read and interpreted by the coding. Agents can also use `PyAutoLens-JAX` database and aggregator tools to combine 
large samples of analysis results into summary tables and figures. This ultimately enables users to analyse and 
interpret modelling results for thousands of lenses entirely through natural language.

# Science projects and open science

When a user begins a specific scientific study using an AI CLI coding agent, `PyAutoLens-Assistant` can create a 
dedicated science project: a logically structured folder linked to a GitHub repository containing the datasets, 
configuration files, analysis scripts, results, plotting scripts and a full transcript with the assistant for 
reproducibility. Every script generated by the assistant is fully documented and can be converted automatically into 
a Jupyter notebook. 
This persistent record of all work allows users to resume work in a later session (which could be with a 
different coding agent) without losing context or progress.
The GitHub repository provides a straightforward way to share results with collaborators, so they can inspect the project’s current state,
understand how each analysis was performed, and provide suggestions or build on the project. Projects can also interface 
with HPC facilities, for example by generating SLURM batch scripts locally for the user to transfer and submit. If the 
study leads to a paper, the completed repository can therefore serve as the paper’s open-source 
companion, enabling readers to reproduce the study end to end or fork it as the starting point for further research.

# Benchmark

We benchmark `PyAutoLens-Assistant` using `ChatGPT` (GPT-5.6 Sol and **GPT-5.5**), `Claude` (Opus 5.0 and **Sonnet 5.0**), 
`Claude Code` (Opus 4.8), `Codex` (GPT-5.6 Sol), `Gemini CLI` (**Gemini 3.1 Pro**), and `OpenCode` (**DeepSeek V4 Flash Free**), 
with bold models available without a paid 
subscription at the time of writing. The primary benchmark uses the COSMOS-Web Ring prompt above, with success 
determined by whether the generated script constructs the requested model and recovers the expected lens configuration. 
All tested models pass, suggesting this task is accessible even to less capable models.

<!--
Further benchmarks reproduce 
the dark-matter subhalo detection in SDSS J0946+1006 [@Vegetti2010; @Minor2021], assess Teacher Mode, and simulate and jointly 
model CCD imaging, interferometric and weak lensing observations of a group-scale lens. The final benchmark intentionally combines 
simulation, multiple data types, a multi-galaxy mass model, and joint inference, which are documented over many different 
`autolens_workspace` examples. Thus it shows how the assistant can generalize across the wikis and skills to do complex
tasks which are not documented individually.
-->

# Model Context Protocol

Model Context Protocol (MCP) provides a standard way for an AI assistant to access external data and software, 
rather than relying only on information copied into a conversation. `PyAutoLens-Assistant` ships with read-only MCP 
tools that allow an AI chat assistant or AI coding agent to search, filter and retrieve completed lens-modelling results 
through natural language. For example, a Euclid lensing analysis may contain models for over 10,000 strong lenses. 
Using MCP, a scientist can point `PyAutoLens-Assistant` to a result hosting server and ask the AI agent 
to “show lenses with Einstein radii above $1.5''$” or “compare the magnification distributions of two samples”. 
They would then receive the relevant parameters and figures without directly interfacing with the `PyAutoLens-JAX` API at all. 
MCP therefore turns a large collection of modelling outputs into an accessible, interactive scientific resource 
that scientists can explore solely via natural language.

# Natural-language development ecosystem

In March 2026, following more than a decade of hand-written software development, `PyAutoLens-JAX` transitioned
to a fully natural-language, agentic-AI development ecosystem called
[`PyAutoScientist`](https://github.com/PyAutoLabs/PyAutoScientist). The ecosystem is organised as a software organism
whose core repositories mirror the roles of human organs:
[`PyAutoBrain`](https://github.com/PyAutoLabs/PyAutoBrain) acts as the reasoning centre, classifying, planning, and
routing tasks through specialist coding agents; [`PyAutoMind`](https://github.com/PyAutoLabs/PyAutoMind) captures
intent by recording plain-English development requirements and tracking them from initial ideas to completed
implementations; and [`PyAutoMemory`](https://github.com/PyAutoLabs/PyAutoMemory) provides long-term scientific memory
through cross-linked literature wikis and verifiable citations. One can therefore conduct software development
entirely through natural language.

# Related software

- `lenstronomy` [@Birrer2018a] is a multi-purpose package for galaxy-scale lens modelling.
- `JAXtronomy` [@Huang2026JAXtronomy] is a `JAX` port of `lenstronomy` supporting automatic differentiation, GPU acceleration, and batched computation.
- `Herculens` [@Galan2022Herculens] provides differentiable, GPU-accelerated lens modelling using `JAX`.
- `GIGA-Lens` [@Gu2022] provides gradient-based, GPU-accelerated Bayesian inference using `TensorFlow` and `JAX`, and now scales across multiple GPU nodes [@Huang2026GIGALens].
- `TinyLensGPU` [@Cao2025] is a `JAX`-based, GPU-accelerated package designed for fast, automated modelling of large galaxy–galaxy strong-lens samples.
- `Gravity.jl` [@Lombardi2024] performs multi-plane lens modelling and Bayesian inference in `Julia`, specialising in point-like and small extended sources.
- `Lenstool` [@Jullo2007] is widely used for cluster-scale lens modelling.
- `GLEE` [@SuyuHalkola2010] is an established package for galaxy- and cluster-scale lens analyses.
- `glafic` [@Oguri2010glafic] provides parametric modelling and simulation tools for galaxy- and cluster-scale gravitational lenses.

# Dependencies

- Core scientific computing: `Python` [@python], `NumPy` [@numpy], `SciPy` [@scipy], `Astropy` [@astropy1; @astropy2], `matplotlib` [@matplotlib], and `Numba` [@numba].
- Scientific utilities: `scikit-image` [@scikit-image], `scikit-learn` [@scikit-learn], `Colossus` [@colossus], `NetworkX` [@networkx], and `SQLite` [@sqlite].
- JAX and GPU computing: `JAX` [@jax], `Optax` [@optax], `Prodigy` [@prodigy], `BlackJAX` [@blackjax], and `TensorFlow Probability` [@tfp].
- Interferometric calculations: `nufftax` [@nufftax], its underlying `FINUFFT` algorithm [@finufft], and the legacy `PyNUFFT` implementation [@pynufft].
- Critical-curve and caustic calculations: `Jax-Zero-Contour` [@jax_zero_contour].
- Non-linear inference: `dynesty` [@dynesty], `emcee` [@emcee], `nautilus` [@nautilus], and `zeus` [@zeus1; @zeus2].
- Posterior analysis and visualization: `corner` [@corner], `anesthetic` [@anesthetic], and `GetDist` [@getdist].
- PyAuto ecosystem components: `PyAutoFit` [@pyautofit] and `PyAutoGalaxy` [@pyautogalaxy].

# AI usage disclosure

Generative AI tools were used to scaffold this manuscript and assist drafting. All scientific claims, 
citations, and prose are reviewed and verified by the authors.

# Acknowledgements

JWN is supported by an STFC/UKRI Ernest Rutherford Fellowship, Project Reference: ST/X003086/1. 
RGH is supported by STFC Opportunities grant ST/T002565/1.

# References
