version 1.0

workflow run_scattered_logistic_ewas {
  input {
    Array[File] beta_matrix_chunks
    File meta_data
    String output_file
    String outcome
    String covariates
    File logistic_ewas_R_script
    String sample_column = "sample_id"
    Int cpu = 4
    Float memory = 8
    Int disk_size = 50
    Int preemptible = 1
    Int maxRetries = 1
  }

  # Run EWAS on each chunk in parallel (scattered)
  scatter (chunk_file in beta_matrix_chunks) {
    call ewas_logistic_chunk {
      input:
        beta_chunk = chunk_file,
        meta_data = meta_data,
        outcome = outcome,
        covariates = covariates,
        logistic_ewas_R_script = logistic_ewas_R_script,
        sample_column = sample_column,
        cpu = cpu,
        memory = memory,
        disk_size = disk_size,
        preemptible = preemptible,
        maxRetries = maxRetries
    }
  }

  output {
    Array[File] ewas_results = ewas_logistic_chunk.chunk_ewas_results
  }

  meta {
    author: "Yash Pershad"
    email: "yash.pershad@vanderbilt.edu"
    description: "Scattered WDL for parallel logistic regression EWAS analysis"
  }
}

task ewas_logistic_chunk {
  input {
    File beta_chunk
    File meta_data
    String outcome
    String covariates
    File logistic_ewas_R_script
    String sample_column
    Int cpu
    Float memory
    Int disk_size
    Int preemptible
    Int maxRetries
  }

  String chunk_name = basename(beta_chunk, ".tsv")
  String output_file = "~{chunk_name}_results.tsv"

  command <<<
    Rscript ~{logistic_ewas_R_script} \
      ~{beta_chunk} \
      ~{meta_data} \
      ~{output_file} \
      ~{outcome} \
      "~{covariates}" \
      ~{sample_column} \
      5000
  >>>

  output {
    File chunk_ewas_results = "~{output_file}"
  }

  runtime {
    docker: "hpoisner/r_phewas:latest"
    memory: memory + " GB"
    disks: "local-disk " + disk_size + " HDD"
    cpu: cpu
    preemptible: preemptible
    maxRetries: maxRetries
  }
}