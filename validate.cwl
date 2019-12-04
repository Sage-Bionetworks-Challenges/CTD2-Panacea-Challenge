#!/usr/bin/env cwl-runner
#
# Validate SC1
#
cwlVersion: v1.0
class: CommandLineTool
baseCommand: [Rscript, /usr/local/bin/validate.R]

hints:
  DockerRequirement:
    dockerPull: docker.synapse.org/syn20968332/scoring_harness:v1

inputs:
  - id: inputfile
    type: File
  - id: template
    type: File

arguments:
  - valueFrom: $(inputs.inputfile.path)
    prefix: --inputfile
  - valueFrom: $(inputs.template.path)
    prefix: --template
  - valueFrom: results.json
    prefix: --results

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
