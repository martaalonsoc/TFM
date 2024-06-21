## Jorge de la Barrera, 2019

workflow VC {

  String sample_id
  String rootdir
  String queue
  File alignment
  File alignment_index

  Object wf = read_json("${rootdir}/0000_pipeline/_json/VC.json")
  String genome_metadata = wf.genome_metadata
  Object genome = read_json("${genome_metadata}")

  #Array[Array[String]] intervals_list = read_tsv("/data_PESA/WGS/ETC/0000_config/intervals/Homo_sapiens_assembly38.intervals")


  call getIntervals {
      input :
        execution = wf.getIntervals,
        queue=queue
    }

  scatter (interval in getIntervals.intervals_array) {
    call VariantCalling {
      input : 
        rootdir = rootdir,
        execution = wf.VariantCalling,
        queue=queue,
        genome = genome,
        sname = sample_id,
        interval = interval,
        alignment = alignment,
        alignment_index = alignment_index
    }
  }

  call MergeVCFs {
    input :
      rootdir = rootdir,
      execution = wf.MergeVCFs,
      queue=queue,
      genome=genome,
      sname = sample_id,
      input_vcfs = VariantCalling.output_gvcf,
      input_vcfs_indexes=VariantCalling.output_gvcf_index
  }
  output {
    File output_gvcf = MergeVCFs.output_gvcf
    File output_gvcf_index = MergeVCFs.output_gvcf_index
    File output_gvcf_md5 = MergeVCFs.output_gvcf_md5
  }
}

task getIntervals {

  Object execution
  String queue
  
  command {
      echo 'test'
  }
  runtime {
    backend : 'SGE_nope'
    cpu : execution.cpu
    memory : execution.memory
    accounting : execution.accounting
    sge_project : queue
  }
  output {
    Array[File] intervals_array = glob("/data_PESA/WGS/ETC/0000_config/scatterIntervals/*.interval_list")
  }
}

task VariantCalling {

  String rootdir
  Object execution
  String queue
  Object genome
  String sname
  String interval
  File alignment
  File alignment_index

  String file_suffix = basename("${interval}", "interval_list")

  command {

    # Load enviroment
    source /programs/GATK/env_vep.sh

    gatk HaplotypeCaller --gatk-config-file ${rootdir}/0000_config/gatk/GATKConfig.properties --java-options "${execution.java_args_HaplotypeCaller}" \
      -R /data_PESA/WGS/REFERENCES/BROAD_hg38/v0/${genome.ref_fasta} \
      -I ${alignment} \
      -O ${sname}.${file_suffix}.vcf.gz \
      -L ${interval} \
      -ERC GVCF \
      -G StandardAnnotation -G AS_StandardAnnotation -G StandardHCAnnotation \
      --max-alternate-alleles 3 \
      --enable-all-annotations true \
      --pair-hmm-implementation ${execution.gatk_gkl_pairhmm_implementation} \
      --native-pair-hmm-threads ${execution.gatk_gkl_pairhmm_threads} \
      --smith-waterman ${execution.smith_waterman_implementation}

      rc=$?;echo "ExitCode:$rc"

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
    File output_gvcf = "${sname}.${file_suffix}.g.vcf.gz"
    File output_gvcf_index = "${sname}.${file_suffix}.g.vcf.gz.tbi"
  }
}

task MergeVCFs {
  String rootdir
  Object execution
  String queue
  Object genome
  String sname
  Array[File] input_vcfs
  Array[File] input_vcfs_indexes

  String outdir = "${rootdir}/_out/VC"

  command {

    # Load enviroment
    source /programs/GATK/env.sh

    mkdir -p "${outdir}"

    picard ${execution.java_args_MergeVCFs} MergeVcfs \
      TMP_DIR=$TMPDIR \
      I=${sep=' I=' input_vcfs} \
      O=${sname}.vcf.gz && \
    gatk ValidateVariants --gatk-config-file ${rootdir}/0000_config/gatk/GATKConfig.properties --java-options "${execution.java_args_ValidateVariants}" \
      -R /data_PESA/WGS/REFERENCES/BROAD_hg38/v0/${genome.ref_fasta} \
      -V ${sname}.vcf.gz \
      -L /data_PESA/WGS/REFERENCES/BROAD_hg38/v0/${genome.wgs_calling_interval_list} \
      --validate-GVCF \
      --create-output-variant-md5 true \
      --dbsnp /data_PESA/WGS/REFERENCES/BROAD_hg38/v0/${genome.dbSNP_vcf} && \
    mv -f ${sname}.vcf.* ${outdir}/

    rc=$?;echo "ExitCode:$rc"

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
    File output_gvcf = "${outdir}/${sname}.g.vcf.gz"
    File output_gvcf_index = "${outdir}/${sname}.g.vcf.gz.tbi"
    File output_gvcf_md5 = "${outdir}/${sname}.g.vcf.gz.md5"
  }
}
