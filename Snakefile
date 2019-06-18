rule all:
    input:
        expand("overview_{type}.pdf", type=['naive', 'corrected']),
        expand("fits_{type}.pdf", type=['naive', 'corrected'])

rule overview_naive:
    input:
        rscript = "overview_naive.r",
        tissues = "tissues.txt",
        orflib = "data/ORF_DMSO_2019-02.txt"
    output:
        outfile = "overview.rds",
        plotfile = "overview_naive.pdf"
    shell:
        "Rscript {input.rscript}"
            " --tissues {input.tissues}"
            " --orflib {input.orflib}"
            " --outfile {output.outfile}"
            " --plotfile {output.plotfile}"

rule overview_corrected:
    input:
        rscript = "overview_corrected.r",
        infile = "overview.rds"
    output:
        plotfile = "overview_corrected.pdf"
    shell:
        "Rscript {input.rscript}"
            " --infile {input.infile}"
            " --plotfile {output.plotfile}"

fields = {
    "naive" : "LFC DMSO/ETP",
    "corrected" : "z_LFC"
}

rule fits:
    input:
        rscript = "fits.r",
        infile = "overview.rds"
    params:
        field = lambda wc : fields[wc.type]
    output:
        outfile = "fits_{type}.xlsx",
        plotfile = "fits_{type}.pdf"
    shell:
        "Rscript {input.rscript}"
            " --infile {input.infile}"
            " --field '{params.field}'"
            " --outfile {output.outfile}"
            " --plotfile {output.plotfile}"
