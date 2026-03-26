version 1.0

# Main workflow for creating methylation matrix from CX_report files
workflow optimized_matrix_generation {
    input {
        String directory
    }

    call find_sorted_files {
        input:
            directory = directory
    }

    call create_methylation_matrix {
        input:
            sorted_methylation_files = find_sorted_files.sorted_files,
    }

    output {
        File beta_matrix = create_methylation_matrix.beta_matrix
    }
}

task find_sorted_files {
    input {
        String directory
        Int cpu = 1
        Float memory = 2
    }

    command <<<
        set -eu -o pipefail
        
        # List all CX_report.txt.gz files
        gsutil ls ~{directory}*sorted.txt > sorted_files.txt
        
        echo "Found $(wc -l < sorted_files.txt) sorted files"
    >>>
    
    output {
        Array[File] sorted_files = read_lines("sorted_files.txt")
    }

    runtime {
        docker: "google/cloud-sdk:latest"
        memory: memory + " GB"
        cpu: cpu
    }
}

