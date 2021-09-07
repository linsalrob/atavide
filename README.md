# atavide

A simple metagenome processing pipeline that is designed to clean the data, assemble it into MAGs, and do something interesting with it.

It is definitely a work in progress, but you can run it with the following command

```bash
snakemake --configfile config/process_metagenomes.yaml -s workflow/process_metagenomes.snakefile --profile slurm
```

But you will need a [slurm profile](https://fame.flinders.edu.au/blog/2021/08/02/snakemake-profiles-updated) to make this work!
