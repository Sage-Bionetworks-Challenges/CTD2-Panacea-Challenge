#!/usr/bin/env cwl-runner
#
# Score SC1
#
cwlVersion: v1.0
class: CommandLineTool
baseCommand: [Rscript, /usr/local/bin/score.R]

hints:
  DockerRequirement:
    dockerPull: docker.synapse.org/syn20968332/scoring_harness:v1

inputs:
  - id: inputfile
    type: File
  - id: goldstandard
    type: File
  - id: nullmodel1
    type: File
  - id: nullmodel2
    type: File
  - id: round
    type: string
  - id: check_validation_finished
    type: boolean?

arguments:
  - valueFrom: $(inputs.inputfile.path)
    prefix: --inputfile
  - valueFrom: $(inputs.goldstandard.path)
    prefix: --goldstandard
  - valueFrom: results.json
    prefix: --results
  - valueFrom: $(inputs.nullmodel1.path)
    prefix: --nullmodel1
  - valueFrom: $(inputs.nullmodel2.path)
    prefix: --nullmodel2
  - valueFrom: $(inputs.round)
    prefix: --round

requirements:
  - class: InlineJavascriptRequirement
     
outputs:
  - id: results
    type: File
    outputBinding:
      glob: results.json
