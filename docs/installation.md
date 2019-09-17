![](../images/ikmb_bfx_logo.png)

# Installation

If you are at the IKMB, you can safely ignore the following...

## References

This pipeline requires the GATK resource bundle as a reference. Please download it from [here](https://software.broadinstitute.org/gatk/download/bundle) to
a location where it is visible to all nodes in your cluster. The folder containing the various assembly-specific sub folders will be needed later. 

## Config file

Several config files are used to provide information on compute setup and locatio of references. 

At the IKMB, you do not have to specify a profile, the default will work. 

Profiles are defined in the file nextflow.config in the `profiles` section. A typical profile will look as follows:

```bash
profiles {
        standard {
                includeConfig 'config/base.config'
                includeConfig 'config/rzcluster.config'
                includeConfig 'config/singularity.config'
                includeConfig 'config/resources.config'
        }
} 
```

The mandatory parts here are `base.config`, which defines the resource requirements of the individual processes, and `resources.config` which specifies
the location of the various resource files needed by the pipelines. 

In addition, one or more config files are needed to define the specs of your cluster and the software provisioning to use. The above example is
set up to use Singularity - replacing the singularity.config with conda.config will use Conda instead. For details on the cluster specific part, have a look at `rzcluster.config`, which is set up for our local SLURM cluster. More information on resource managers, please see the Nextflow [documentation](https://www.nextflow.io/docs/latest/executor.html).

Adding all of this to the config file may then look like so:

```bash

profiles {
        standard {
                includeConfig 'config/base.config'
                includeConfig 'config/rzcluster.config'
                includeConfig 'config/singularity.config'
                includeConfig 'config/resources.config'
        }
        your_cluster {
                includeConfig 'config/base.config'
                includeConfig 'config/your_cluster.config'
                includeConfig 'config/singularity.config'
                includeConfig 'config/resources.config'
        }
}

```

### Resource bundle

Most importantly, the pipeline relies on the GATK resource bundle (see above) - the location of the folder containing the various assemblies
must be specified in your cluster config file with `params.gatk_bundle_path`.


