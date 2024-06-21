## Jorge de la Barrera, 2021

workflow BQSR {
#BQSR bins the qualities which makes a significantly smaller bam
  String sample_id=sample_id,
  String queue=queue,
  String accounting=accounting,
  Object wf=wf.BQSR,
  File genome_metadata
  File sequence_grouping_file
  File sequence_grouping_with_unmapped_file
  File alignment
  File alignment_index

  Object genome = read_json(genome_metadata)
  Array[Array[String]] sequence_grouping = read_tsv(sequence_grouping_file)
  Array[Array[String]] sequence_grouping_with_unmapped = read_tsv(sequence_grouping_with_unmapped_file)

  scatter (subgroup in sequence_grouping) {
    call BaseRecalibrator {
      input :
        rootdir = rootdir,
        execution = wf.BaseRecalibrator,
        queue=queue,
        genome = genome,
        sname = sample_id,
        sequence_group_interval = subgroup,
        alignment = alignment,
        alignment_index = alignment_index
    }
  }
  call GatherBqsrReports {
    input : 
      rootdir = rootdir,
      execution = wf.GatherBqsrReports,
      queue=queue,
      genome = genome,
      sname = sample_id,
      input_bqsr_reports = BaseRecalibrator.table
  }

  scatter (subgroup in sequence_grouping_with_unmapped) {
    call ApplyBQSR {
      input :
        rootdir = rootdir,
        execution = wf.ApplyBQSR,
        queue=queue,
        genome = genome,
        sname = sample_id,
        recal_data = GatherBqsrReports.bqsr_report,
        sequence_group_interval = subgroup,
        alignment = alignment,
        alignment_index = alignment_index
    }
  }
  call GatherBamFiles {
    input : 
      rootdir = rootdir,
      execution = wf.GatherBamFiles,
      queue=queue,
      sname = sample_id,
      alignment_in = ApplyBQSR.out_alignment
  }
  output {
    File out_alignment = GatherBamFiles.alignment
    File out_alignment_index = GatherBamFiles.alignment_index
  }
}


task BaseRecalibrator {

  String rootdir
  Object execution
  String queue
  Object genome
  String sname
  Array[String] sequence_group_interval
  File alignment
  File alignment_index

  Array[String] known_sites = prefix("/data_PESA/WGS/REFERENCES/BROAD_hg38/v0/", genome.known_indels_sites_VCFs)

  command {

    # Load enviroment
    source /programs/GATK/env.sh

    # Debug
    sh ${rootdir}/0000_pipeline/_debug.sh

    gatk BaseRecalibrator --gatk-config-file ${rootdir}/0000_config/gatk/GATKConfig.properties --java-options "${execution.java_args_BaseRecalibrator}" \
    -R /data_PESA/WGS/REFERENCES/BROAD_hg38/v0/${genome.ref_fasta} \
    -I ${alignment} \
    --use-original-qualities \
    -O ${sname}.recal_data.table \
    --known-sites /data_PESA/WGS/REFERENCES/BROAD_hg38/v0/${genome.dbSNP_vcf} \
    --known-sites ${sep=' --known-sites ' known_sites} \
    -L ${sep=" -L " sequence_group_interval}

    rc=$?

    # Saca el RC
    echo "ExitCode:$rc"

    exit $rc
  }
  runtime {
    backend : 'SGE_nope'
    cpu : execution.cpu
    memory : execution.memory
    accounting : execution.accounting
    sge_project : queue
  }
  output {
    File table = "${sname}.recal_data.table"
  }
}

task GatherBqsrReports {

  String rootdir
  Object execution
  String queue
  Object genome
  String sname
  Array[File] input_bqsr_reports

  command {

    # Load enviroment
    source /programs/GATK/env.sh

    # Debug
    sh ${rootdir}/0000_pipeline/_debug.sh

    gatk GatherBQSRReports --gatk-config-file ${rootdir}/0000_config/gatk/GATKConfig.properties --java-options "${execution.java_args_GatherBqsrReports}" \
        -I ${sep=' -I ' input_bqsr_reports} \
        -O ${sname}.recal_data.csv

    rc=$?

    # Saca el RC
    echo "ExitCode:$rc"

    exit $rc
  }
  runtime {
    backend : 'SGE_nope'
    cpu : execution.cpu
    memory : execution.memory
    accounting : execution.accounting
    sge_project : queue
  }
  output {
    File bqsr_report = "${sname}.recal_data.csv"
  }
}

task ApplyBQSR {

  String rootdir
  Object execution
  String queue
  Object genome
  String sname
  File recal_data
  Array[String] sequence_group_interval
  File alignment
  File alignment_index

  Array[String] known_sites = prefix("/data_PESA/WGS/REFERENCES/BROAD_hg38/v0/", genome.known_indels_sites_VCFs)

  command {

    # Load enviroment
    source /programs/GATK/env.sh

    # Debug
    sh ${rootdir}/0000_pipeline/_debug.sh

    gatk ApplyBQSR --gatk-config-file ${rootdir}/0000_config/gatk/GATKConfig.properties --java-options "${execution.java_args_ApplyBQSR}" \
      --add-output-sam-program-record \
      -R /data_PESA/WGS/REFERENCES/BROAD_hg38/v0/${genome.ref_fasta} \
      -I ${alignment} \
      --use-original-qualities \
      -O ${sname}.bam \
      --create-output-bam-index \
      -bqsr ${recal_data} \
      --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 \
      -L ${sep=" -L " sequence_group_interval} && \
    samtools quickcheck -q ${sname}.bam

    rc=$?

    # Saca el RC
    echo "ExitCode:$rc"

    exit $rc
  }
  runtime {
    backend : 'SGE_nope'
    cpu : execution.cpu
    memory : execution.memory
    accounting : execution.accounting
    sge_project : queue
  }
  output {
    File out_alignment = "${sname}.bam"
  }
}

task GatherBamFiles {

  String rootdir
  Object execution
  String queue
  String sname
  Array[File] alignment_in
  
  String indir = "${rootdir}/_out/sortBAM"
  String outdir = "${rootdir}/_out/BQSR"

  command {

    # Load enviroment
    source /programs/GATK/env.sh
   
    # Debug
    sh ${rootdir}/0000_pipeline/_debug.sh

    #Copy results
    mkdir -p "${outdir}"

    picard ${execution.java_args_GatherBamFiles} GatherBamFiles \
      TMP_DIR=$TMPDIR \
      INPUT=${sep=' INPUT=' alignment_in} \
      OUTPUT=${sname}.bam \
      CREATE_INDEX=true \
      CREATE_MD5_FILE=false && \
    samtools quickcheck -q ${sname}.bam && \
    mv -f ${sname}.ba* ${outdir}/
    
    rc=$?

    # Saca el RC
    echo "ExitCode:$rc"

    if [ $rc -eq '0' ]; then 
      rm -rf ../../call-BaseRecalibrator && \
      rm -rf ../../call-GatherBqsrReports && \
      rm -rf ../../call-ApplyBQSR && \
      rm ${indir}/${sname}.ba* && \
      rmdir ${indir}
    fi

    exit $rc
  }
  runtime {
    backend : 'SGE_nope'
    cpu : execution.cpu
    memory : execution.memory
    accounting : execution.accounting
    sge_project : queue
  }
  output {
    File alignment = "${outdir}/${sname}.bam"
    File alignment_index = "${outdir}/${sname}.bai"
  }
}