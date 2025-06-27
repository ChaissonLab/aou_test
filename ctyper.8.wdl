version 1.0
workflow RunCtyperBatch {
  input {
    Array[File] input_bams
    Array[File] input_bai_files
    Array[String] output_tab_filenames
    File matrix
    File matrixIndex
    File matrixBin
    File ref
    File fai
  }

  call CtyperTask {
    input:
      ref=ref,
      fai=fai,
      bam_files = input_bams,
      bam_bai_files = input_bai_files,
      output_tab_filenames = output_tab_filenames,
      matrix = matrix,
      matrixIndex = matrixIndex,
      matrixBin = matrixBin,
      n_threads=8
  }

  output {
    Array[File] tabs = CtyperTask.output_files
  }
}

task CtyperTask {
  input {
    Array[File] bam_files
    Array[File] bam_bai_files
    Array[String] output_tab_filenames
    File matrix
    File matrixIndex
    File matrixBin
    File ref
    File fai
    Int n_threads
  }
  output {
   Array[File] output_files = output_tab_filenames
  }

  command <<<
    # Write list of expected VCF output names to a file
    /bin/ctyper2 -T ~{ref} -I ~{write_lines(bam_files)} -O ~{write_lines(output_tab_filenames)} -m ~{matrix} -n 8
  >>>


  runtime {
    docker: "mchaisso/ctyper2:v3"
    memory: "64G"
    cpu: 8
    preemptible: 3
    disks: "local-disk 1000 HDD"
  }
}

