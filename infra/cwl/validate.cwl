#!/usr/bin/env cwl-runner
#
# Validate SC1
#
cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - Rscript
  - /usr/local/bin/validation_scoring.R

hints:
  DockerRequirement:
    dockerPull: docker.synapse.org/syn20968332/TO_UPDATE_WITH_TAG

inputs:
  - id: inputfile
    type: File
  - id: goldstandard
    type: File
  - id: entity_type
    type: string

arguments:
  - valueFrom: $(inputs.inputfile)
    prefix: -s
  - valueFrom: $(inputs.goldstandard.path)
    prefix: -g
  - valueFrom: $(inputs.entity_type)
    prefix: -e
  - valueFrom: results.json
    prefix: -r

requirements:
  - class: InlineJavascriptRequirement
     
outputs:
  - id: results
    type: File
    outputBinding:
      glob: results.json   

  - id: status
    type: string
    outputBinding:
      glob: results.json
      loadContents: true
      outputEval: $(JSON.parse(self[0].contents)['prediction_file_status'])

  - id: invalid_reasons
    type: string
    outputBinding:
      glob: results.json
      loadContents: true
      outputEval: $(JSON.parse(self[0].contents)['prediction_file_errors'])
