3#!/bin/bash

# Run deformetrica
source activate deformetrica
deformetrica estimate model.xml data_set.xml --p optimization_parameters.xml -v DEBUG
conda deactivate deformetrica
