# Some key papers in species distribution modelling

This paper by Stephen Phillips talks about the problem with failing to account for spatial sampling bias. I also introduces the target group background method: <https://doi.org/10.1890/07-2153.1>

This paper by Corey Merow gives a practical guide to using Maxent. It talks about the FactorBiasOut method (bias layer method) that we used: <https://doi.org/10.1111/j.1600-0587.2013.07872.x>


This paper by Will Fithian has another method for accounting for bias, that uses a small set of presence-absence data in conjunction to correct presence-only data for spatial sampling bias:
<https://doi.org/10.1111/2041-210X.12242>

This paper by Guru Guillera-Arroita talks about when an SDM is ‘good enough’ for the purpose you are using it for - what you can and can’t do with different kinds of SDMs:
<https://doi.org/10.1111/geb.12268>

This paper is beyond the scope of our course; it’s a bit technical and less practical so don’t worry about reading it unless you are interested. It explains why Poisson point processes are the natural way of modelling presence-only data (and why MaxEnt is actually fitting a Poisson point process model): <https://doi.org/10.1111/2041-210X.12352>

And this other paper by Guru explains why you shouldn’t ever use maxent for modelling presence-absence data! (apparently some people do this): <https://doi.org/10.1111/2041-210X.12252>

This paper discusses the different types of validation metrics that are available for presence-only and presence-absence data, and which are appropriate for which purpose: <https://doi.org/10.1111/2041-210X.12123>

This paper by Jane Elith is a very practical and user-friendly introduction to using Boosted Regression Trees (a nice machine learning method) with a working application to species distribution modelling
<https://doi.org/10.1111/j.1365-2656.2008.01390.x>

This paper by Roozbeh Valavi compares how well a range of presence-only methods predict, and shows how you go about setting up and fitting your can make a difference
<https://doi.org/10.1002/ecm.1486>

This other paper of Jane Elith’s talks about how to model species that are in the process of spreading through a landscape: <https://doi.org/10.1111/j.2041-210X.2010.00036.x>

And this *other* paper by Jane Elith describes the statistical model of maxent and talks about how  decisions you need to make in fitting your model in maxent. It's a must read for using maxent for SDMs.
<https://doi.org/10.1111/j.1472-4642.2010.00725.x>

This paper by many leading lights in the SDM space describes a standardising protocol for reporting species distribution analyses --- the ODMAP protocol. <https://doi.org/10.1111/ecog.04960>
