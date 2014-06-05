humanomics
=======

Pipelines for analyses of -omics experiments in human

Usage
-------

For now, just run:

./0--ConfigurePipeline.bash

./0--MakePipeline.py

This will start by writing a configuration file (pipeline.conf, by default) which will be used to produce the following SLURM scripts:

```0--FastQC.bash```

1--MappingAndPreProcessing.bash

2a--QualityControl.bash

2b--HaplotypeCaller.bash

