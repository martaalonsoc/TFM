## Jorge de la Barrera, 2019

workflow sortBAM {

  String sample_id
  String rootdir
  String queue
  File alignment

  Object wf = read_json("${rootdir}/0000_pipeline/_json/sortBAM.json")
  Object execution = wf.execution

  call mappedBAM_2_sortedBAM {
    input :
      rootdir = rootdir,
      execution = execution,
      queue=queue,
      sname = sample_id,
      in_alignment = alignment
  }

  output {
    File out_alignment = mappedBAM_2_sortedBAM.alignment
    File out_alignment_index = mappedBAM_2_sortedBAM.alignment
  }
}

task mappedBAM_2_sortedBAM {

  String rootdir
  Object execution
  String queue
  String sname
  File in_alignment
  
  String indir = "${rootdir}/_out/mapping"
  String outdir = "${rootdir}/_out/sortBAM"

  command {

    # Load enviroment
    source /programs/GATK/env.sh

    ulimit -n 50000
   
    # Debug
    sh ${rootdir}/0000_pipeline/_debug.sh

    mkdir -p "${outdir}"

    # job
    samtools sort \
        -T $TMPDIR \
        -m ${execution.samtools_memory_by_thread} \
        --threads ${execution.cpu} \
        -l ${execution.compression_level} \
        ${in_alignment} -o ${sname}.bam && \
    samtools index \
        -@ ${execution.cpu} \
        ${sname}.bam ${sname}.bai && \
    samtools quickcheck -q ${sname}.bam && \
    mv -f ${sname}.ba* ${outdir}/

    rc=$?

    # Saca el RC
    echo "ExitCode:$rc"    

    if [ $rc -eq '0' ]; then 
      rm ${indir}/${sname}.ba* && \
      rmdir ${indir}
    fi

    exit $rc   

}
  runtime {
    backend : 'SGE_pthreads'
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
