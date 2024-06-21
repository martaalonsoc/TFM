## Jorge de la Barrera, 2019

workflow BAM2CRAM {

  String sample_id
  String rootdir
  String queue
  File alignment
  File alignment_index

  Object wf = read_json("${rootdir}/0000_pipeline/BAM2CRAM.json")
  String genome_metadata = wf.genome_metadata
  Object genome = read_json(genome_metadata)
  Object execution = wf.execution

  call toCRAM {
    input :
      rootdir = rootdir,
      execution = execution,
      queue=queue,
      genome = genome,
      sname = sample_id,
      in_alignment = alignment,
      in_alignment_index = alignment_index
  }

  output {
    File out_alignment = mappedBAM_2_sortedBAM.alignment
    File out_alignment_index = mappedBAM_2_sortedBAM.alignment
  }
}

task toCRAM {

  String rootdir
  Object execution
  String queue
  Object genome
  String sname
  File in_alignment
  File in_alignment_index
 
  String basedir = "${rootdir}/${sname}"
  String outdir = "${basedir}/_out/BAM2CRAM"

  command {

    # Load enviroment
    source /programs/GATK/env.sh

    # Debug
    sh ${rootdir}/0000_pipeline/_debug.sh

    # job
    samtools view --threads ${execution.cpu} -C \
    -T /data_PESA/WGS/REFERENCES/BROAD_hg38/v0/${genome.ref_fasta} \
      ${in_alignment} > ${sname}.cram && \
      samtools quickcheck -q ${sname}.cram

    rc=$?

    # Loga si esta montado el home (a stdout)
    echo -ne "automount_home:`mount -l|grep "/home/jdelabarrera"`" > /dev/stdout

    # Saca el RC
    echo "ExitCode:$rc"    

    # Copy results
    mkdir -p "${outdir}"
    ln -f ${sname}.bam ${outdir}/${sname}.cram
    ln -f ${sname}.bai ${outdir}/${sname}.bai

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
