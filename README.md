# popsims
## A simulation on population genetics

>[!NOTE]
>This simulation is designed for practical sessions in Population Genetics courses (MBG211 & MBG451).

This agent-based simulation uses Mesa module.

## Running the simulation ###

To initiate the simulation simply execute run_popsim.py. Enter the values for initial parameters (See below.).

E.g. you can use `%run run_popsim.py` on Jupyter Notebook.

The prompt will request a *key* to initiate the simulation. The *visual* key initiates the visual simulation and the *run* key initiates simulation and plots frequencies over generations.

>[!WARNING]
> 09.01.25 : The visual simulation is currently not working properly. Please use "run" key.

 `Enter the frequency of dominant allele (p) :` asks for the frequency of wild-type allele, *p* , to estimate the mutant allele frequency, *q*, and then the initial genotypic frequencies.

 E.g. p = 0.7

 `Enter the mutability (probably for random mutation into recessive allele) : ` asks for the probability of a random mutation creating the mutant allele.
 E.g. mutability = 0.01

 `Enter reproduction rate :` asks for the ratio of population that reproduces and creates an offspring every generation. It should be between 0 and 1.
 E.g. reproduction rate = 0.2

 `Enter selection coefficient for homozygous dominant genotype :`
 
 `Enter selection coefficient for heterozygous genotype :`

 `Enter selection coefficient for homozygous recessive genotype :`

 These prompts asks for selection coefficient for the genotypes. E.g. if there is selection against the recessive phenotype with s=0.4, then you can enter 0 for first, 0 for second and 0.4 for third prompt.
 
`Enter initial population size : ` asks for the initial population size. E.g. N = 1000.

`Enter width :` and `Enter height :` creates a map of *width* x *height* squares with an organism on each square. Therefore you should enter width and height such that their product is higher then the initial population size. E.g. 40 and 40.

`Enter number of generations :` asks for the duration of the simulation. E.g. 100




## Initial parameters 
**p**                  : The frequency of the dominant allele. q and genotypic frequences of initial population is estimated over p.

**mut_prob**           : Probability for a dominant WT allele mutating into the recessive mutant allele.

**survival_prob**      : Probability of survival (1-selection_coefficient) for each genotype.

**pop_size**           : Initial population size

**grid**               : Creates the map of the habitat.

**track_frequencies**  : Tracks changes in genotypic and allelic frequencies each generation.

**n_alive**            : Number of alive organisms each generation.
## Version 0 -- Simple Hardy Weinberg Equilibrium with Natural Selection and Mutations

Each **organism** in the population is introduced as an *agent* and the **population** dynamics is introduced as the *model*.

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



**r_move**             : Radius of available locations for an offspring. Introduces "resource limit".

**N_gen**              : Generations for repeated runs.
