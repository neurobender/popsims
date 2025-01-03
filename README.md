# popsims
## A simulation on population genetics

[!NOTE]
This simulation is designed for practical sessions in Population Genetics courses (MBG211 & MBG451).

This agent-based simulation uses Mesa module.

Each **organism** in the population is introduced as an *agent* and the **population** dynamics is introduced as the *model*.

## Version 0 -- Simple Hardy Weinberg Equilibrium with Natural Selection and Mutations

The population dynamics are defined with a expanded *Hardy-Weinberg Equilibrium* model with following assumptions :
* A single gene with two alleles, complete/incomplete/co-dominance.
  * Wild-type (WT) (**+**) is dominant over the mutant (**m**).
  * Genotypes are defined as homozygous WT/dominant (**++**), heterozygous (**+m**) and homozygous recessive/mutant (**mm**). 
* Random mating.
  * The probability of reproduction is capped with a **reproduction rate** , representing the ratio of reproducing individuals in a generation over the total population size.
* No migration.
* No population subdivision.

**Random genetic drift** is introduced as the accumulated impact of the following factors :
* Randomized individual reproduction event.
* Random mate selection.
  * Each organism randomly chooses a mate of any available genotypes.
* Random allele segregation and fertilization.
  * The gamete that will be used in fertilization may contain any of the two alleles from each parent.
* Randomized survival with limiting conditions if there is natural selection.
  * Selection coefficients for each genotype is externally set.

**Mutability** is introduced for the prevalence of a random mutation from dominant WT --> recessive mutant allele.

**Crowdedness** is introduced as a radius of available area with resources around a parent.

## Initial parameters 
p                  : The frequency of the dominant allele. q and genotypic frequences of initial population is estimated over p.

mut_prob           : Probability for a dominant WT allele mutating into the recessive mutant allele.

survival_prob      : PRobability of survival (1-selection_coefficient) for each genotype.

pop_size           : Initial population size

grid               : Creates the map of the habitat.

track_frequencies  : Tracks changes in genotypic and allelic frequencies each generation.

n_alive            : Number of alive organisms each generation.

r_move             : Radius of available locations for an offspring. Introduces "resource limit".

N_gen              : Generations for repeated runs.
