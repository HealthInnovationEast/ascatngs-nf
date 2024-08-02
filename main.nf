#!/usr/bin/env nextflow
nextflow.enable.dsl=2

def helpMessage() {
    // TODO
    log.info """
    Please see here for usage information: https://github.com/HealthInnovationEast/ascatngs-nf/blob/main/docs/Usage.md
    """.stripIndent()
}

// Show help message
if (params.help) {
  helpMessage()
  exit 0
}

/*--------------------------------------------------------
  Defining and showing header with all params information
----------------------------------------------------------*/

// Header log info

def summary = [:]

if (workflow.revision) summary['Pipeline Release'] = workflow.revision

summary['Output dir']                                  = params.outdir
summary['Launch dir']                                  = workflow.launchDir
summary['Working dir']                                 = workflow.workDir
summary['Script dir']                                  = workflow.projectDir
summary['User']                                        = workflow.userName

// then arguments
summary['pairs']                                       = params.pairs
summary['genome_fa']                                   = params.genome_fa
summary['gender_snps']                                 = params.gender_snps
summary['snp_gc_corr']                                 = params.snp_gc_corr

log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "-\033[2m--------------------------------------------------\033[0m-"

// Importantly, in order to successfully introspect:
// - This needs to be done first `main.nf`, before any (non-head) nodes are launched.
// - All variables to be put into channels in order for them to be available later in `main.nf`.

ch_repository         = Channel.of(workflow.manifest.homePage)
ch_commitId           = Channel.of(workflow.commitId ?: "Not available is this execution mode. Please run 'nextflow run ${workflow.manifest.homePage} [...]' instead of 'nextflow run main.nf [...]'")
ch_revision           = Channel.of(workflow.manifest.version)
ch_scriptName         = Channel.of(workflow.scriptName)
ch_scriptFile         = Channel.of(workflow.scriptFile)
ch_projectDir         = Channel.of(workflow.projectDir)
ch_launchDir          = Channel.of(workflow.launchDir)
ch_workDir            = Channel.of(workflow.workDir)
ch_userName           = Channel.of(workflow.userName)
ch_commandLine        = Channel.of(workflow.commandLine)
ch_configFiles        = Channel.of(workflow.configFiles)
ch_profile            = Channel.of(workflow.profile)
ch_container          = Channel.of(workflow.container)
ch_containerEngine    = Channel.of(workflow.containerEngine)

/*----------------------------------------------------------------
  Setting up additional variables used for documentation purposes
-------------------------------------------------------------------*/

Channel
    .of(params.raci_owner)
    .set { ch_raci_owner }

Channel
    .of(params.domain_keywords)
    .set { ch_domain_keywords }

/*----------------------
  Setting up input data
-------------------------*/

// Define Channels from input
// only if not in dsl2

/*-----------
  Processes
--------------*/

process ascat_counts {
    input:
        path('genome.fa')
        path('genome.fa.fai')
        path('snp.gc')
        path('sex.loci')
        tuple val(groupId), val(type), val(sampleId), val(protocol), val(platform), file(htsfile), file(htsidx), val(purity), val(ploidy)

    output:
        tuple path("${sampleId}.count.gz"), path("${sampleId}.count.gz.tbi")
        path("${sampleId}.is_male.txt")
        tuple val(groupId), val(type), val(sampleId), val(protocol), val(platform), path("${sampleId}.count.gz"), path("${sampleId}.count.gz.tbi"), val(purity), val(ploidy), path("${sampleId}.is_male.txt"), emit: to_ascat

    // makes sure pipelines fail properly, plus errors and undef values
    shell = ['/bin/bash', '-euo', 'pipefail']

    stub:
        """
        echo ">>>ascat_counts htsfile: ${htsfile}"
        touch ${sampleId}.count.gz
        touch ${sampleId}.count.gz.tbi
        touch ${sampleId}.is_male.txt
        """

    script:
        """
        # remove logs for sucessful jobs
        export PCAP_THREADED_REM_LOGS=1
        ascatCounts.pl -o . \
            -b $htsfile \
            -r genome.fa \
            -sg snp.gc \
            -l sex.loci \
            -c $task.cpus
        """
}

process ascat {
    input:
        path('genome.fa')
        path('genome.fa.fai')
        path('genome.fa.dict')
        path('snp.gc')
        tuple val(groupId), val(types), val(sampleIds), val(protocol), val(platform), path(counts), path(indexes), val(purity), val(ploidy), path(isMale)

    output:
        tuple path('*.copynumber.caveman.vcf.gz'), path('*.copynumber.caveman.vcf.gz.tbi')
        path('*.png')
        tuple val(groupId), path('*.copynumber.caveman.csv'), path('*.samplestatistics.txt'), emit: ascat_for_caveman
        path('*.copynumber.txt.gz')
        tuple path("*.count.gz", includeInputs: true), path("*.count.gz.tbi", includeInputs: true)
        path('*.is_male.txt', includeInputs: true)

    publishDir {
        def case_idx = types.indexOf('case')
        def ctrl_idx = types.indexOf('control')
        "${params.outdir}/${sampleIds[case_idx]}_vs_${sampleIds[ctrl_idx]}"
    }, mode: 'copy'

    // makes sure pipelines fail properly, plus errors and undef values
    shell = ['/bin/bash', '-euo', 'pipefail']

    stub:
        def case_idx = types.indexOf('case')
        def ctrl_idx = types.indexOf('control')
        def purity_val = purity[case_idx] != 'NA' ? "--purity ${purity[case_idx]}" : ''
        def ploidy_val = ploidy[case_idx] != 'NA' ? "--ploidy ${ploidy[case_idx]}" : ''

        """
        echo ">>>ascat purity: ${purity_val}"
        echo ">>>ascat ploidy: ${ploidy_val}"
        echo ">>>ascat tumour: ${counts[case_idx]}"
        echo ">>>ascat normal: ${counts[ctrl_idx]}"
        echo ">>>ascat tumour: ${isMale[case_idx]}"
        echo ">>>ascat normal: ${isMale[ctrl_idx]}"
        touch ${sampleIds[case_idx]}_vs_${sampleIds[ctrl_idx]}.copynumber.caveman.vcf.gz
        touch ${sampleIds[case_idx]}_vs_${sampleIds[ctrl_idx]}.copynumber.caveman.vcf.gz.tbi
        touch ${sampleIds[case_idx]}_vs_${sampleIds[ctrl_idx]}.png
        touch ${sampleIds[case_idx]}_vs_${sampleIds[ctrl_idx]}.copynumber.caveman.csv
        touch ${sampleIds[case_idx]}_vs_${sampleIds[ctrl_idx]}.samplestatistics.txt
        touch ${sampleIds[case_idx]}_vs_${sampleIds[ctrl_idx]}.copynumber.txt.gz
        # others are from inputs
        """

    script:
        def case_idx = types.indexOf('case')
        def ctrl_idx = types.indexOf('control')
        def purity_val = purity[case_idx] != 'NA' ? "--purity ${purity[case_idx]}" : ''
        def ploidy_val = ploidy[case_idx] != 'NA' ? "--ploidy ${ploidy[case_idx]}" : ''

        """
        # remove logs for sucessful jobs
        export PCAP_THREADED_REM_LOGS=1
        SPECIES=`head -n 2 genome.fa.dict | tail -n 1 | perl -ne 'm/SP:([^\t]+)/;print \$1;'`
        ASSEMBLY=`head -n 2 genome.fa.dict | tail -n 1 | perl -ne 'm/AS:([^\t]+)/;print \$1;'`
        ascat.pl -nb -f -o . \
            ${purity_val} \
            ${ploidy_val} \
            -ra "\$ASSEMBLY" -rs "\$SPECIES" \
            -pr "${protocol[ctrl_idx]}" -pl "${platform[ctrl_idx]}" \
            -r genome.fa \
            -sg snp.gc \
            -g ${isMale[ctrl_idx]} \
            -t ${counts[case_idx]} -tn ${sampleIds[case_idx]} \
            -n ${counts[ctrl_idx]} -nn ${sampleIds[ctrl_idx]} \
            -c $task.cpus
        """
}

workflow {
    genome_fa    = file(params.genome_fa)
    genome_fai    = file("${params.genome_fa}.fai")
    genome_dict    = file("${params.genome_fa}.dict")
    gender_snps  = file(params.gender_snps)
    snp_gc_corr  = file(params.snp_gc_corr)

    if(params.counts) {
        // input is different to allow for count files and is_male files
        pairs = Channel.fromPath(params.pairs)
        case_control_map = pairs.splitCsv(header: true).map {
            row -> tuple(row.groupId, row.type, row.sampleId, row.protocol, row.platform, file(row.reads), file(row.readIdx), row.purity, row.ploidy, file(row.isMale))
        }
    }
    else {
        pairs = Channel.fromPath(params.pairs)
        case_control_map = pairs.splitCsv(header: true).map {
            row -> tuple(row.groupId, row.type, row.sampleId, row.protocol, row.platform, file(row.reads), file(row.readIdx), row.purity, row.ploidy)
        }
    }

    main:

        if(params.counts) {
            to_ascat = case_control_map
        }
        else {
            ascat_counts(
                genome_fa,
                genome_fai,
                snp_gc_corr,
                gender_snps,
                case_control_map
            )
            to_ascat = ascat_counts.out.to_ascat
        }

        ascat(
            genome_fa,
            genome_fai,
            genome_dict,
            snp_gc_corr,
            to_ascat.groupTuple()
        )
}
