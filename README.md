This repo contains a few scripts we use to perform optimization for dadi (https://bitbucket.org/gutenkunstlab/dadi). Current contents:

2popIM.py: fits a 2-population joint SFS to a model with the following parameters:
	a population split time
	the size of population one after the split
	the size of population two after the split
	population one's present-day size (exponential change immediately after split)
	population two's present-day size (exponential change immediately after split)
	the rate of migration from population two into population one
	the rate of migration from population one into population two
	ancestral theta is also extracted and used for scaling
	usage: python 2popIM.py inputSFSFile swarmSize generationsPerYear
		inputSFSFile is a two-population SFS properly formatted for dadi
		swarmSize is the number of CPUs to use for particle swarm optimization
		generationsPerYear is used for scaling the split time from generations into years

2popIsolation.py: same as 2popIM.py but with no migration

2popIM\_sym\_fixed\_ancGrowth.py: similar to 2popIM.py, but with symmetric migration, and also an epoch of exponential change in ancestral pop
	this model also contains a fixed rate of SNP mispolarization (currently 0.01).

2popISO\_pMisFixed\_ancGrowth.py: same as 2popIM\_sym\_fixed\_ancGrowth.py but with no migration

dadiFunctions.py: contains various demographic functions that can be used for dadi's integration

Note that pyOpt is also required to use these scripts (http://www.pyopt.org/index.html).
