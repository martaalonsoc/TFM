## Jorge de la Barrera, 2019

workflow mapping {

  String sample_id
  String queue
  String accounting
  Object wf
  File genome_metadata
  Array[File] reads_in

  Object genome = read_json(genome_metadata)

  scatter (sample_reads_in in reads_in) {
    call uBAM_2_mappedBAM {
      input :
        execution = wf.uBAM_2_mappedBAM,
        queue = queue,
        accounting = accounting,
        ref_fasta = genome.basedir + genome.bwa_fasta,
        reads_in = sample_reads_in
    }
  }

  call MarkDuplicates {
    input :
      sample_id =sample_id,
      execution = wf.MarkDuplicates,
      queue = queue,
      accounting = accounting,
      alignment_in =  uBAM_2_mappedBAM.alignment
  }

  call mappedBAM_2_sortedBAM {
    input :
      sample_id =sample_id,
      execution = wf.MarkDuplicates,
      queue = queue,
      accounting = accounting,
      alignment_in =  MarkDuplicates.alignment
  }

  output {
    Array[File] bwa_stderr = uBAM_2_mappedBAM.bwa_stderr
    File markduplicates_metrics = MarkDuplicates.markduplicates_metrics
    File out_alignment = mappedBAM_2_sortedBAM.alignment
    File out_alignment_index = mappedBAM_2_sortedBAM.alignment_index
  }
}


task uBAM_2_mappedBAM {

  Object execution
  String queue
  String accounting
  String ref_fasta
  String reads_in

  String file_prefix = basename(reads_in, "_unmapped.bam")
  String bwa_commandline="-K 100000000 -p -Y -v 3 -t ${execution.cpu} ${ref_fasta}"

  command {

    # Load enviroment
    source /programs/GATK/env.sh

    # job
    picard ${execution.java_args_SamToFastq} SamToFastq \
      TMP_DIR=$TMPDIR \
      I=${indir}/${file_prefix}.unmapped.bam \
      F=/dev/stdout \
      INTERLEAVE=true \
      NON_PF=true | \
    bwa mem ${bwa_commandline} /dev/stdin - 2>${file_prefix}.bwa.stderr.log | \
    picard ${execution.java_args_MergeBamAlignment} MergeBamAlignment \
        TMP_DIR=$TMPDIR \
        VALIDATION_STRINGENCY=SILENT \
        ALIGNED_BAM=/dev/stdin \
        UNMAPPED_BAM=${indir}/${file_prefix}.unmapped.bam \
        R=${ref_fasta} \
        O=${file_prefix}.bam \
        EXPECTED_ORIENTATIONS=FR \
        ATTRIBUTES_TO_RETAIN=X0 \
        ATTRIBUTES_TO_REMOVE=NM \
        ATTRIBUTES_TO_REMOVE=MD \
        ATTRIBUTES_TO_REMOVE=UQ \
        PAIRED_RUN=true \
        SORT_ORDER="unsorted" \
        IS_BISULFITE_SEQUENCE=false \
        ALIGNED_READS_ONLY=false \
        CLIP_ADAPTERS=false \
        MAX_RECORDS_IN_RAM=2000000 \
        ADD_MATE_CIGAR=true \
        MAX_INSERTIONS_OR_DELETIONS=-1 \
        PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
        UNMAPPED_READ_STRATEGY=COPY_TO_TAG \
        ALIGNER_PROPER_PAIR_FLAGS=true \
        UNMAP_CONTAMINANT_READS=true \
        ADD_PG_TAG_TO_READS=false && \   
    samtools quickcheck -q ${file_prefix}.bam && \
    grep -m1 "read .* ALT contigs" ${file_prefix}.bwa.stderr.log | grep -v "read 0 ALT contigs" 
  }
  runtime {
    backend : 'SGE_nope'
    cpu : execution.cpu
    memory : execution.memory
    accounting : accounting
    sge_project : queue
  }
  output {
    File alignment = "${file_prefix}.bam"
  }
}

task MarkDuplicates {

  String sample_id
  Object execution
  String queue
  String accounting
  Array[File] alignment_in
  Int? optical_dup_pixel_dist=2500
  # El valor por defecto es valido para casi todos los casos. Si se indica null se desactiva 
  # la deteccion de duplicados opticos
  String? read_name_regex

  command {

    # Load enviroment
    source /programs/GATK/env.sh

    picard ${execution.java_args_MarkDuplicates} MarkDuplicates \
      TMP_DIR=$TMPDIR \
      I=${sep=' INPUT=' alignment_in} \
      O=${sname}.bam \
      M=${sname}.markduplicates_metrics \
      VALIDATION_STRINGENCY=SILENT \
      ${" --READ_NAME_REGEX=" + read_name_regex} \
      OPTICAL_DUPLICATE_PIXEL_DISTANCE=${optical_dup_pixel_dist} \
      ASSUME_SORT_ORDER="queryname" \
      CLEAR_DT=false \
      ADD_PG_TAG_TO_READS=false \
      MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=20000 && \
    samtools quickcheck -q ${sname}.bam
    rc=$?

    # Saca el RC
    echo "ExitCode:$rc"

    if [ $rc -eq '0' ]; then 
      rm -rf ../../call-uBAM_2_mappedBAM 
    fi
        
    exit $rc
  }
  runtime {
    backend : 'SGE_nope'
    cpu : execution.cpu
    memory : execution.memory
    accounting : accounting
    sge_project : queue
  }
  output {
    File alignment = "${outdir}/${sname}.bam"
  }
}

task mappedBAM_2_sortedBAM {

  String sample_id
  Object execution
  String queue
  String accounting
  File alignment_in

  command {

    # Load enviroment
    source /programs/GATK/env.sh

    # job
    samtools sort \
        -T $TMPDIR \
        -m ${execution.samtools_memory_by_thread} \
        --threads ${execution.cpu} \
        -l ${execution.compression_level} \
        ${alignment_in} -o ${sample_id}.bam && \
    samtools index \
        -@ ${execution.cpu} \
        ${sample_id}.bam ${sample_id}.bai && \
    samtools quickcheck -q ${sample_id}.bam
    rc=$?

    # Saca el RC
    echo "ExitCode:$rc"    

    exit $rc   

}
  runtime {
    backend : 'SGE_pthreads'
    cpu : execution.cpu
    memory : execution.memory
    accounting : accounting
    sge_project : queue
  }
  output {
    File alignment = "${outdir}/${sname}.bam"
    File alignment_index = "${outdir}/${sname}.bai"
  }
}


