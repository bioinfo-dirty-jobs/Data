from subprocess import run

if len(snakemake.input) > 2:
    cmd = [snakemake.params['cmd'],
           '-p',
           '-T', str(snakemake.threads),
           '-t', snakemake.params['gtf_feature_type'],
           '-a', snakemake.params['gtf'][0],
           '-g', snakemake.params['gtf_attribute_type'],
           '-o', snakemake.output[0],
           snakemake.input['bam']]
else:
    cmd = [snakemake.params['cmd'],
           '-T', str(snakemake.threads),
           '-t', snakemake.params['gtf_feature_type'],
           '-a', snakemake.params['gtf'],
           '-g', snakemake.params['gtf_attribute_type'],
           '-o', snakemake.output[0],
           snakemake.input['bam']]
run(cmd)
