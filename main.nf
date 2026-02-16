nextflow.enable.dsl=2

process HELLO {
  output:
    stdout

  """
  echo "Hello from Nextflow"
  """
}

workflow { HELLO() }
