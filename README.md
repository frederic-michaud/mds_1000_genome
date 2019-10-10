# mds_1000_genome

This repository contain a snakemake script to download the 1000 genome project, filter the data and run medeas on them. 

It can be run using: 
`snakemake --directory local_data --cores 80`

where the configuration files have been moved to the local_data directory

A representation of the data can be seen using:

`snakemake --dag | dot -Tpdf > dag.pdf`

A conda environment in which the pipeline works can be created using:

`conda env create -n mds_1000_genome -f conda_envs.yaml`

`source activate mds_1000_genome`

It can be run using the cluster.json on a slurm cluster using:

`Snakemake -j 80 --directory local_data --cluster-config cluster.json --cluster "sbatch -p {cluster.partition} --mem={cluster.mem} -t {cluster.time} -A {cluster.account} --cpus-per-task={cluster.cpus} --job-name={cluster.name} --error={cluster.error} --output={cluster.output}"`
