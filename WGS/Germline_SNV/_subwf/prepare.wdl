## Jorge de la Barrera, 2019

workflow prepare {

  String sample_id
  String rootdir
  String queue

  Object wf = read_json("${rootdir}/0000_pipeline/_json/prepare.json")
  
  Object sample = read_json("/data_PESA/WGS/ETC/0000_dataset/all_samples/${sample_id}.json")

  String sample_name = sample.name
  Array[Object] files = sample.files

  scatter (file in files) {

    Array[String] fastq = file.fastq

    call GetRG_hashCode {
     input :
          execution = wf.GetRG_hashCode,
          queue = queue,
          fastq = fastq
    }

    call Fastq_2_uBAM {
      input :
        rootdir = rootdir,
        execution = wf.Fastq_2_uBAM,
        queue = queue,
        hashCode = GetRG_hashCode.hashCode,
        sname = sample_name,
        sdata = file
    }
  }
}

# Obtiene @RG hash
task GetRG_hashCode {

  Object execution
  String queue
  Array[String] fastq

  String file_prefix = basename(fastq[0], ".fastq.gz")

  command {
        echo "${file_prefix}"|sha1sum -b|head -c 6
  }
  runtime {
    backend : 'SGE_nope'
    cpu : execution.cpu
    memory : execution.memory
    accounting : execution.accounting
    sge_project : queue
  }
  output {
    String hashCode = read_string(stdout())
  }
}

task Fastq_2_uBAM {

  String rootdir
  Object execution
  String queue
  String hashCode
  String sname
  Object sdata

  Array[String] fastq = sdata.fastq
  String file_prefix = basename(fastq[0], ".fastq.gz")
  
  String logdir = "${rootdir}/_log/prepare"
  String outdir = "${rootdir}/_out/prepare"

  command {

    # Load enviroment
    source /programs/GATK/env.sh

    # Debug
    sh ${rootdir}/0000_pipeline/_debug.sh
  
    # Create output directory 
    mkdir -p "${logdir}"
    mkdir -p "${outdir}"

    picard ${execution.java_args_FastqToSam} FastqToSam \
      TMP_DIR=$TMPDIR \
      F1=/data_PESA/WGS/RAW/${sdata.fastq_path}/${fastq[0]} \
      F2=/data_PESA/WGS/RAW/${sdata.fastq_path}/${fastq[1]} \
      O=/dev/stdout \
      SM=${sname} \
      RG=${sname}.${hashCode} \
      PU=${sdata.flowcell}.${sdata.lane}.${sdata.barcode} \
      DT=${sdata.date} \
      LB=${sdata.library} \
      PL=${sdata.platform_technology} \
      PM=${sdata.platform_model} \
      CN=${sdata.sequencing_center} | \
      picard ${execution.java_args_MarkIlluminaAdapters} MarkIlluminaAdapters \
        TMP_DIR=$TMPDIR \
        I=/dev/stdin \
        O=${file_prefix}.unmapped.bam \
        M=${file_prefix}_markilluminaadapters_metrics.txt && \
      mv ${file_prefix}.unmapped.bam ${outdir}/${file_prefix}.unmapped.bam


      rc=$?

      # Saca el RC
      echo "ExitCode:$rc"

      # Copy results
      
      mv ${file_prefix}_markilluminaadapters_metrics.txt ${logdir}/${file_prefix}_markilluminaadapters_metrics.txt

      exit $rc
  }
  runtime {
    backend : 'SGE_nope'
    cpu : execution.cpu
    memory : execution.memory
    accounting : execution.accounting
    sge_project : queue
  }
}

