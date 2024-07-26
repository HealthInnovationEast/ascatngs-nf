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
summary['core_ref']                                    = params.core_ref
summary['cvn_sv']                                      = params.cvn_sv
summary['qc_genotype']                                 = params.qc_genotype

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

// Do not delete this process
// Create introspection report

process obtain_pipeline_metadata {
    publishDir "${params.tracedir}", mode: "copy"

    input:
      val(repository)
      val(commit)
      val(revision)
      val(script_name)
      val(script_file)
      val(project_dir)
      val(launch_dir)
      val(work_dir)
      val(user_name)
      val(command_line)
      val(config_files)
      val(profile)
      val(container)
      val(container_engine)
      val(raci_owner)
      val(domain_keywords)

    output:
      path("pipeline_metadata_report.tsv"), emit: pipeline_metadata_report

    // same as script except ! instead of $ for variables
    shell:
      '''
      echo "Repository\t!{repository}"                  > temp_report.tsv
      echo "Commit\t!{commit}"                         >> temp_report.tsv
      echo "Revision\t!{revision}"                     >> temp_report.tsv
      echo "Script name\t!{script_name}"               >> temp_report.tsv
      echo "Script file\t!{script_file}"               >> temp_report.tsv
      echo "Project directory\t!{project_dir}"         >> temp_report.tsv
      echo "Launch directory\t!{launch_dir}"           >> temp_report.tsv
      echo "Work directory\t!{work_dir}"               >> temp_report.tsv
      echo "User name\t!{user_name}"                   >> temp_report.tsv
      echo "Command line\t!{command_line}"             >> temp_report.tsv
      echo "Configuration file(s)\t!{config_files}"    >> temp_report.tsv
      echo "Profile\t!{profile}"                       >> temp_report.tsv
      echo "Container\t!{container}"                   >> temp_report.tsv
      echo "Container engine\t!{container_engine}"     >> temp_report.tsv
      echo "RACI owner\t!{raci_owner}"                 >> temp_report.tsv
      echo "Domain keywords\t!{domain_keywords}"       >> temp_report.tsv
      awk 'BEGIN{print "Metadata_variable\tValue"}{print}' OFS="\t" temp_report.tsv > pipeline_metadata_report.tsv
      '''

    stub:
      '''
      touch pipeline_metadata_report.tsv
      '''
}

process prep_ref {
    input:
        file(core_ref)
        file(cvn_sv)
        file(qc_genotype)

    output:
        path 'ref', type: 'dir', emit: ref
        path 'ascat/SnpGcCorrections.tsv', emit: snps_gc
        path 'gender.tsv', emit: snps_sex
        path'ref_cache', type: 'dir', emit: ref_cache

    shell = ['/bin/bash', '-euo', 'pipefail']

    stub:
        """
        mkdir -p ref
        touch ref/genome.{fa,fa.fai,fa.dict}
        mkdir -p ascat
        touch ascat/SnpGcCorrections.tsv
        touch gender.tsv
        mkdir -p ref_cache
        """

    script:
        """
        mkdir ref
        tar --strip-components 1 -C ref -zxvf $core_ref
        tar --strip-components 1 --no-anchored -zxvf $cvn_sv SnpGcCorrections.tsv
        tar --strip-components 1 --no-anchored -zxvf $qc_genotype gender.tsv

        # build ref-cache
        seq_cache_populate.pl -subdirs 2 -root ./ref_cache ref/genome.fa
        """
}

process ascat_counts {
    input:
        path('ref')
        path('snp.gc')
        path('sex.loci')
        tuple val(groupId), val(type), val(sampleId), val(protocol), val(platform), file(htsfile), file(htsidx), file(htsStats)

    output:
        tuple path("${sampleId}.count.gz"), path("${sampleId}.count.gz.tbi")
        path("${sampleId}.is_male.txt")
        tuple val(groupId), val(type), val(sampleId), val(protocol), val(platform), path("${sampleId}.count.gz"), path("${sampleId}.count.gz.tbi"), path("${sampleId}.is_male.txt"), emit: to_ascat

    // makes sure pipelines fail properly, plus errors and undef values
    shell = ['/bin/bash', '-euo', 'pipefail']

    stub:
        """
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
            -r ref/genome.fa \
            -sg snp.gc \
            -l sex.loci \
            -c $task.cpus
        """
}

process ascat {
    input:
        path('ref')
        path('snp.gc')
        tuple val(groupId), val(types), val(sampleIds), val(protocol), val(platform), path(counts), path(indexes), path(ismale)

    output:
        tuple path('*.copynumber.caveman.vcf.gz'), path('*.copynumber.caveman.vcf.gz.tbi')
        path('*.png')
        tuple val(groupId), path('*.copynumber.caveman.csv'), path('*.samplestatistics.txt'), emit: ascat_for_caveman
        path('*.copynumber.txt.gz')
        tuple path("*.count.gz", includeInputs: true), path("*.count.gz.tbi", includeInputs: true)

    publishDir {
        def case_idx = types.indexOf('case')
        def ctrl_idx = types.indexOf('control')
        "${params.outdir}/${sampleIds[case_idx]}_vs_${sampleIds[ctrl_idx]}/ascat"
    }, mode: 'copy'

    // makes sure pipelines fail properly, plus errors and undef values
    shell = ['/bin/bash', '-euo', 'pipefail']

    stub:
        def case_idx = types.indexOf('case')
        def ctrl_idx = types.indexOf('control')
        """
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
        """
        # remove logs for sucessful jobs
        export PCAP_THREADED_REM_LOGS=1
        SPECIES=`head -n 2 ref/genome.fa.dict | tail -n 1 | perl -ne 'm/SP:([^\t]+)/;print \$1;'`
        ASSEMBLY=`head -n 2 ref/genome.fa.dict | tail -n 1 | perl -ne 'm/AS:([^\t]+)/;print \$1;'`
        ascat.pl -nb -f -o . \
            -ra "\$ASSEMBLY" -rs "\$SPECIES" \
            -pr "${protocol[ctrl_idx]}" -pl "${platform[ctrl_idx]}" \
            -r ref/genome.fa \
            -sg snp.gc \
            -g ${ismale[ctrl_idx]} \
            -t ${counts[case_idx]} -tn ${sampleIds[case_idx]} \
            -n ${counts[ctrl_idx]} -nn ${sampleIds[ctrl_idx]} \
            -c $task.cpus
        """
}

workflow {
    core_ref     = file(params.core_ref)
    cvn_sv       = file(params.cvn_sv)
    qc_genotype  = file(params.qc_genotype)

    pairs = Channel.fromPath(params.pairs)

    case_control_map = pairs.splitCsv(header: true).map { row -> tuple(row.groupId, row.type, row.sampleId, row.protocol, row.platform, file(row.reads), file(row.readIdx), file(row.readStats)) }

    main:
        obtain_pipeline_metadata(
            ch_repository,
            ch_commitId,
            ch_revision,
            ch_scriptName,
            ch_scriptFile,
            ch_projectDir,
            ch_launchDir,
            ch_workDir,
            ch_userName,
            ch_commandLine,
            ch_configFiles,
            ch_profile,
            ch_container,
            ch_containerEngine,
            ch_raci_owner,
            ch_domain_keywords
        )

        prep_ref(
            core_ref,
            cvn_sv,
            qc_genotype
        )

        ascat_counts(
            prep_ref.out.ref,
            prep_ref.out.snps_gc,
            prep_ref.out.snps_sex,
            case_control_map
        )
        ascat(
            prep_ref.out.ref,
            prep_ref.out.snps_gc,
            ascat_counts.out.to_ascat.groupTuple()
        )

    /*
    STUB run
    rm -rf results/* .nextflow* work
    nextflow run main.nf \
        -profile test -stub-run \
        --core_ref data/core_ref_GRCh38_hla_decoy_ebv.tar.gz \
        --cvn_sv data/CNV_SV_ref_GRCh38_hla_decoy_ebv_brass6+.tar.gz \
        --qc_genotype data/qcGenotype_GRCh38_hla_decoy_ebv.tar.gz \
        --pairs data/test.csv

    rm -rf results/* .nextflow* work
    nextflow run $HOME/git/HealthInnovationEast/ascatngs-nf/main.nf \
        -profile test,singularity,slurm -resume \
        --core_ref data/core_ref_GRCh38_hla_decoy_ebv.tar.gz \
        --cvn_sv data/CNV_SV_ref_GRCh38_hla_decoy_ebv_brass6+.tar.gz \
        --qc_genotype data/qcGenotype_GRCh38_hla_decoy_ebv.tar.gz \
        --pairs /home/kr525/git/cynapse-ccri/cgpwgs-nf/test.csv
    */

}
