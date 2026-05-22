#!/bin/bash
#SBATCH -c 11                              # Request cores
#SBATCH -t 0-12:00                         # Runtime in D-HH:MM format
#SBATCH --partition gpu_quad               # Use a quad GPU
#SBATCH --gres=gpu:rtx8000:1,vram:26G      # Number to use
#SBATCH --mem-per-cpu=10G                   # Memory per core
#SBATCH -o jobs/deform_fafb_%j.out           # File to which STDOUT will be written, including job ID (%j)
#SBATCH -e jobs/deform_fafb_%j.err           # File to which STDERR will be written, including job ID (%j)

start=`date +%s`

source /home/ab714/bancpipeline/o2/o2_env.sh
ulimit -c 0  # suppress core dumps from R shutdown crashes

cd /n/data1/hms/neurobio/wilson/banc/deformetrica/fafb/Neuron/
source activate deformetrica
deformetrica estimate model.xml data_set.xml --p optimization_parameters.xml -v DEBUG || true
conda deactivate

cd /n/data1/hms/neurobio/wilson/banc/deformetrica/fafb/NeuronMesh/
source activate deformetrica
deformetrica estimate model.xml data_set.xml --p optimization_parameters.xml -v DEBUG || true
conda deactivate

cd /n/data1/hms/neurobio/wilson/banc/deformetrica/fafb/Points/
source activate deformetrica
deformetrica estimate model.xml data_set.xml --p optimization_parameters.xml -v DEBUG || true
conda deactivate

cd /n/data1/hms/neurobio/wilson/banc/deformetrica/fafb/Surface/
source activate deformetrica
deformetrica estimate model.xml data_set.xml --p optimization_parameters.xml -v DEBUG || true
conda deactivate

cd /n/data1/hms/neurobio/wilson/banc/deformetrica/fafb/
source activate deformetrica
deformetrica estimate model.xml data_set.xml --p optimization_parameters.xml -v DEBUG || true
conda deactivate

Rscript $HOME/bancpipeline/deform/transfer_vtk.R

end=`date +%s`
runtime=$((end-start))
echo "script completed in: "
echo $runtime