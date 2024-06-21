## Jorge de la Barrera, 2021

import "mapping.wdl" as mapping
import "BQSR.wdl" as BQSR
import "FixValidate.wdl" as FixValidate
import "VC.wdl" as VC

workflow Germline_SNV_sample {

	String sample_id
  String queue
  String accounting
  Object wf
  File genome_metadata
  Array[File] reads_in


  call mapping.mapping as mapping {
    input : 
      sample_id=sample_id,
      queue=queue,
      accounting=accounting,
      wf=wf.mapping,
      genome_metadata=genome_metadata,
      reads_in=reads_in
  }

  call BQSR.BQSR as BQSR {
    input : 
      sample_id=sample_id,
      queue=queue,
      accounting=accounting,
      wf=wf.BQSR,
      genome_metadata=genome_metadata,
      sequence_grouping_file=wf.sequence_grouping_file.
      sequence_grouping_with_unmapped_file=wf.sequence_grouping_with_unmapped_file,
      alignment= mapping.out_alignment,
      alignment_index=mapping.out_alignment_index     
  }

  call FixValidate.FixValidate as FixValidate {
    input : 
      sample_id=sample_id,
      rootdir=rootdir,
      queue=queue,
      alignment= BQSR.out_alignment,
      alignment_index= BQSR.out_alignment_index
  }

  call VC.VC as VC {
    input : 
      sample_id=sample_id,
      rootdir=rootdir,
      queue=queue,
      alignment= FixValidate.out_alignment,
      alignment_index= FixValidate.out_alignment_index
  }

  output {
    #mapping
    Array[File] bwa_stderr = mapping.bwa_stderr
    File markduplicates_metrics = mapping.markduplicates_metrics
    # Alineamiento
    File alignment = FixValidate.out_alignment
    File alignment_index = FixValidate.out_alignment_index
    File alignment_md5 = FixValidate.out_alignment_md5
    File validation_report = FixValidate.validation_report

    # gVCF
    File output_gvcf = VC.output_gvcf
    File output_gvcf_index = VC.output_gvcf_index
    File output_gvcf_md5 = VC.output_gvcf_md5
  }
}