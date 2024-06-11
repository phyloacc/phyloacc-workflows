VALUES = [i/10.0 for i in range(1,10)]

OUTDIR = "/n/holylfs05/LABS/informatics/Users/gthomas/phyloacc-workflows/scripts/hmm-benchmarks/"

rule all:
    input: 
        expand("conserved-elements-20-{t0_0}-{t1_1}-{e0_0}-{e1_1}-1.0.bed", t0_0=VALUES, t1_1=VALUES, e0_0=VALUES, e1_1=VALUES)

rule hmm_predict:        
    output: 
        bed = "conserved-elements-20-{t0_0}-{t1_1}-{e0_0}-{e1_1}-1.0.bed"
    params:
        outdir = OUTDIR
    log:
        log = "conserved-elements-20-{t0_0}-{t1_1}-{e0_0}-{e1_1}-1.0.log"
    resources:
        partition = "shared",
        time = "00:35:00",
        mem = "1G"
    shell:
        """
        cd {params.outdir}
        python ../conserved_elements_hmm.py {wildcards.t0_0} {wildcards.t1_1} {wildcards.e0_0} {wildcards.e1_1} 1.0
        """
