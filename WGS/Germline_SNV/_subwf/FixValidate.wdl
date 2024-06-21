## Jorge de la Barrera, 2019

workflow FixValidate {

  String sample_id
  String rootdir # Aqui indica el directorio flowcell/muestra ('sampledir')
  String queue
  File alignment
  File alignment_index
  
  String sample_name = basename(alignment, ".bam")

  Object wf = read_json("${rootdir}/0000_pipeline/_json/FixValidate.json")
  String genome_metadata = wf.genome_metadata
  Object genome = read_json(genome_metadata)

  call SetNmAndUqTags {
    input :
      rootdir = rootdir,
      execution = wf.SetNmAndUqTags,
      queue=queue,
      sname = sample_name,
      alignment_in = alignment,
      alignment_index_in = alignment_index,
      genome = genome
  }

  call ValidateSAM {
    input :
      rootdir = rootdir,
      execution = wf.ValidateSAM,
      queue=queue,
      sname = sample_name,
      alignment_in =  SetNmAndUqTags.alignment,
      alignment_index_in =  SetNmAndUqTags.alignment_index,
      alignment_md5_in =  SetNmAndUqTags.alignment_md5,
      genome = genome
  }
  output {
    File out_alignment = ValidateSAM.alignment
    File out_alignment_index = ValidateSAM.alignment_index
    File out_alignment_md5 = ValidateSAM.alignment_md5
    File validation_report = ValidateSAM.validation_report
  }
}


task SetNmAndUqTags {

  String rootdir
  Object execution
  String queue
  String sname
  File alignment_in
  File alignment_index_in
  Object genome
 
  command {

    # Load enviroment
    source /programs/GATK/env.sh

    # Debug
    sh ${rootdir}/0000_pipeline/_debug.sh

    # job
    picard ${execution.java_args_SetNmAndUqTags} SetNmAndUqTags \
      TMP_DIR=$TMPDIR \
      I=${alignment_in} \
      O=${sname}.bam \
      CREATE_INDEX=true \
      CREATE_MD5_FILE=true \
      COMPRESSION_LEVEL=5 \
      REFERENCE_SEQUENCE=/data_PESA/WGS/REFERENCES/BROAD_hg38/v0/${genome.ref_fasta} && \
    samtools quickcheck -q ${sname}.bam

  # Saca el RC
    echo "ExitCode:$rc"

    rc=$?

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
    File alignment = "${sname}.bam"
    File alignment_index = "${sname}.bai"
    File alignment_md5 = "${sname}.bam.md5"
  }
}

task ValidateSAM {

  String rootdir
  Object execution
  String queue
  String sname
  File alignment_in
  File alignment_index_in
  File alignment_md5_in
  Object genome

  String indir = "${rootdir}/_out/BQSR"
  String outdir = "${rootdir}/_out/FixValidate" 

  command {

    # Load enviroment
    source /programs/GATK/env.sh
    
    # Debug
    sh ${rootdir}/0000_pipeline/_debug.sh

    mkdir -p "${outdir}"

    picard ${execution.java_args_ValidateSamFile} ValidateSamFile \
      TMP_DIR=$TMPDIR \
      MAX_OPEN_TEMP_FILES=60000 \
      I=${alignment_in} \
      OUTPUT= ${sname}.report \
      REFERENCE_SEQUENCE=/data_PESA/WGS/REFERENCES/BROAD_hg38/v0/${genome.ref_fasta} \
      MODE=VERBOSE \
      IS_BISULFITE_SEQUENCED=false && \
    grep -q "^No errors found$" ${sname}.report && \
    cp -Lf ${alignment_in} ${outdir}/ && \
    cp -Lf ${alignment_index_in} ${outdir}/ && \
    cp -Lf ${alignment_md5_in} ${outdir}/ && \
    mv -f ${sname}.report ${outdir}/ && \
    rm -rf ../../call-SetNmAndUqTags && \
    rm ${indir}/${sname}.ba* && \
    rmdir ${indir}
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
    File alignment_md5 = "${outdir}/${sname}.bam.md5"
    File validation_report = "${outdir}/${sname}.report"
  }
}


