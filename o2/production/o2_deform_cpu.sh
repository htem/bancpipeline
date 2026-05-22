#!/bin/bash
#SBATCH -c 1                              # Request cores
#SBATCH -t 0-12:00                        # Runtime in D-HH:MM format
#SBATCH -p short                          # Partition to run in
#SBATCH --mem-per-cpu=100G                # Memory per core
#SBATCH -o jobs/deform_cpu%j.outb         # File to which STDOUT will be written, including job ID (%j)
#SBATCH -e jobs/deform_cpu%j.err          # File to which STDERR will be written, including job ID (%j)

start=`date +%s`

source /home/ab714/bancpipeline/o2/o2_env.sh
ulimit -c 0  # suppress core dumps from R shutdown crashes

cd /n/data1/hms/neurobio/wilson/banc/deformetrica/fafb/Neuron/
source activate deformetrica
deformetrica estimate model.xml data_set.xml --p optimization_parameters.xml --output='cpu_output/' -v DEBUG || true
conda deactivate

cd /n/data1/hms/neurobio/wilson/banc/deformetrica/fafb/NeuronMesh/
source activate deformetrica
deformetrica estimate model.xml data_set.xml --p optimization_parameters.xml --output='cpu_output/' -v DEBUG || true
conda deactivate

cd /n/data1/hms/neurobio/wilson/banc/deformetrica/fafb/Points/
source activate deformetrica
deformetrica estimate model.xml data_set.xml --p optimization_parameters.xml --output='cpu_output/' -v DEBUG || true
conda deactivate

cd /n/data1/hms/neurobio/wilson/banc/deformetrica/fafb/Surface/
source activate deformetrica
deformetrica estimate model.xml data_set.xml --p optimization_parameters.xml --output='cpu_output/' -v DEBUG || true
conda deactivate

cd /n/data1/hms/neurobio/wilson/banc/deformetrica/fafb/
source activate deformetrica
deformetrica estimate model.xml data_set.xml --p optimization_parameters.xml --output='cpu_output/' -v DEBUG || true
conda deactivate

Rscript $HOME/bancpipeline/deform/transfer_vtk.R

end=`date +%s`
runtime=$((end-start))
echo "script completed in: "
echo $runtime