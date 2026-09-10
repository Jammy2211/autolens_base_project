<!-- Archived 2026-09-03 from README.md (section "#### AI Chat Assistant") (autolens_assistant#120). Marked UNSUPPORTED on
     2026-09-10: conversational chat routes were retired under the agentic-only support policy and are not
     coming back as supported onboarding. Kept for reference only; links inside may point at other archived pages. -->

> **Archived and unsupported (2026-09-10).** Historical text; see the [README](../../README.md) for the supported
> agentic setup.

#### AI Chat Assistant

Depending on which assistant you chose in the guide, your first message may need to open with its setup
instructions (for example the GitHub sync bootstrap prompt). 

With that in place, here is a good initial prompt to try it out, noting that data for the COSMOS-Web Ring is included in this repository as an example:

```
[Setup instructions for your chosen assistant, if its guide page says you need them]

Find the data on the Cosmos-Web ring, give me a short script to plot it in PyAutoLens and then given that I'm a 
new user give me an overview of the different ways we can perform strong lens modeling of this system.
```

The above prompt will give an overview of the PyAutoLens API for plotting, describe how you can perform lens modeling 
of the system, and ask you follow up questions which will get a discussion going so you can begin using PyAutoLens
for a more specific task.

The `autolens_assistant` can easily handle more complex tasks: using the prompt below you'll get an end-to-end Python script for multi-wavelength lens modeling of the COSMOS-Web Ring!

```
[Setup instructions for your chosen assistant, if its guide page says you need them]

I want to model the F277W and F444W JWST imaging of the COSMOS-Web Ring independently, which are in 
the folder dataset/imaging/cosmos_web_ring. Model the lens light with a multi-Gaussian expansion (MGE), its mass with a singular 
isothermal ellipsoid plus external shear, and model the source also using an MGE. For speed, run the analysis on my 
laptop GPU using a JAX optimizer that estimates only the maximum-likelihood solution. Plot the observed image at 
each wavelength in the left column, its lensed source model in the middle column, and its source on the right column.
```

