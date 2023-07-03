# Bayesian-MSC

Bayesian implementation of the Multispecies Coalescent model (MSC) for two species.

The MSC tracks the genealogy of a sample of sequences back through time and, in this way, models how an allele may have arisen from a common ancestor. Two types of parameters are present in the MSC: speciation times ($\tau's$) and population size parameters ($\theta's$). Here, $\theta = 4N\mu$ is the mutational distance separating two randomly sampled sequences in a population, with N referring to the effective population size and $\mu$ the mutation rate per site per generation. For example, $\theta$ = 0.001 means that two sequences will have, on average, ~1 difference per kilobase. Both $\tau$ and $\theta$ are measured using this mutational distance. As we have only one sequence per population, the $\theta$ parameters for the modern populations are unidentifiable. 
