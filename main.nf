#!/usr/bin/env nextflow
nextflow.enable.dsl=2

def helpMessage() {
    // TODO
    log.info """
    Please see here for usage information: https://github.com/cynapse-ccri/template-nf/blob/main/docs/Usage.md
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

// then arguments ******* EDIT THIS
summary['hello']                                       = params.hello


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

process trivial_example {
    input:
        val(hello_thing)

    output:
        path 'greeting.txt', emit: greeting

    shell = ['/bin/bash', '-euo', 'pipefail']

    // this just creates files required to allow testing of a structure
    stub:
        """
        touch greeting.txt
        """

    // this is the work that does what you'd expect
    script:
        """
        echo "Hello ${hello_thing}" > greeting.txt
        """
}

workflow {
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
        trivial_example(params.hello)
}
