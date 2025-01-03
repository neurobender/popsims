# popsims
A simulation on population genetics

This simulation is designed for practical sessions in Population Genetics courses (MBG211 & MBG451).

The agent-based simulation uses Mesa module.

Each organism in the population is introduced as an agent and the population dynamics is introduced as the model.

v0 -- Modified Hardy-Weinberg Population Genetics model with following assumptions :
** A single gene with two alleles, complete dominance. Wild-type is dominant over the mutant.
** Random mating 
** No migration
** No population subdivision

Random genetic drift is introduced with
-- randomized reproduction event
-- random mate selection
-- random allele segregation
-- random fertilization
-- randomized survival with limiting conditions if there is natural selection.

Mutability is introduced for the prevalence of a random mutation from dominant WT --> recessive mutant

Reproduction rate is introduced as a limiting condition.

Crowdedness is introduced as a radius of available area with resources around a parent.

Initial parameters -- external 
p                  : The frequency of the dominant allele. q and genotypic frequences of initial population is estimated over p.
mut_prob           : Probability for a dominant WT allele mutating into the recessive mutant allele.
survival_prob      : PRobability of survival (1-selection_coefficient) for each genotype.
pop_size           : Initial population size
grid               : Creates the map of the habitat.
track_frequencies  : Tracks changes in genotypic and allelic frequencies each generation.
n_alive            : Number of alive organisms each generation.
r_move             : Radius of available locations for an offspring. Introduces "resource limit".
N_gen              : Generations for repeated runs.
