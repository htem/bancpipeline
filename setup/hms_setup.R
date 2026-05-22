# srun --pty -p interactive -t 0-6:00 --mem 25G -c 1 bash

# # Load modules
# module purge
# module load gcc/14.2.0
# module load R/4.4.2
# module load cmake/3.31.2
# module load java/jdk-23.0.1
# # Export udunits variables
# export UDUNITS2_INCLUDE=/n/app/udunits/2.2.28-gcc-9.2.0/include
# export UDUNITS2_LIBS=/n/app/udunits/2.2.28-gcc-9.2.0/lib
# # Start R
# R

### Start R session
# screen -S R
# srun --pty -p interactive -t 0-10:00 --mem 5G -c 1 bash
# module load gcc/9.2.0 R/4.3.1
# R

## RStudio
# https://o2portal.rc.hms.harvard.edu 

# Where data is stored
# /n/data1/hms/neurobio/wilson/fafbz

## Set up python environment
# srun --pty -p interactive -t 0-10:00 --mem 5G -c 1 bash
# module load gcc/6.2.0 miniconda3/4.10.3 python/3.7.4
# conda update conda
# conda  create --name flywire_env python=3.7
# conda init bash
# source /n/app/miniconda3/4.10.3/etc/profile.d/conda.sh
# conda activate flywire_env
# conda install nomkl
# conda install -c conda-forge pykdtree
# pip install meshparty cloud-volume~=1.20.1
# pip3 install git+git://github.com/schlegelp/skeletor@master
# pip3 install fastremap
# pip3 install ncollpyde
# Could not install packages due to an EnvironmentError: [Errno 13] Permission denied: 'f2py3.7'
# Consider using the `--user` option or check the permissions.

## Get PCRE2
# wget https://github.com/PCRE2Project/pcre2/releases/download/pcre2-10.42/pcre2-10.42.tar.bz2
# tar xvf pcre2-10.42.tar.bz2
# cd pcre2-10.42
# ./configure --enable-utf8 --prefix=$HOME
# make
# make install

### Install our own version of R
# wget http://cran.rstudio.com/src/base/R-4/R-4.3.1.tar.gz
# tar xvf R-4.3.1.tar.gz
# cd R-4.3.1
# ./configure --prefix=$HOME/R --with-x=no
# make && make install

## Add these items to our PATH
#touch ~/.bash_profile
#nano ~/.bash_profile
#### Add this line: export PATH="$PATH:/home/hy92/.local/share/r-miniconda/:/home/hy92/.local/share/r-miniconda/bin/:/home/hy92/.local/share/r-miniconda/envs:/home/hy92/.local/share/r-miniconda/condabin:/home/hy92/ninja/include:/home/hy92/fuse-3.16.1/lib:/home/hy92/.local/usr/local/lib/:/home/hy92/rclone-v1.63.1-linux-amd64/:/home/hy92/R-4.3.1/bin/:/home/hy92/.local/bin/:/home/hy92/.local/lib/:/home/hy92/.local/include/:/home/hy92/blender-3.6.1-linux64/:/home/hy92/R/x86_64-pc-linux-gnu-library/:/home/hy92/jdk-11.0.8/bin/:/home/hy92/.R/gargle/gargle-oauth/:/home/hy92/rclone-v1.54.1-linux-amd64/:/home/hy92/opt/libgit2-1.1.0/:/home/hy92/pcre2/pcre2-10.42/:/home/hy92/pcre2/bin/:/home/hy92/elastix-5.1.0/bin:/home/hy92/google-cloud-sdk/bin"
# export LD_LIBRARY_PATH ="$LD_LIBRARY_PATH:/home/hy92/elastix-5.1.0/lib"
# FILES="/files/Neurobio/wilsonlab/"
# DATA="/n/data1/hms/neurobio/wilson"
# TP="/files/Neurobio/wilsonlab/abates/bergIII/virmen_cone"
#### Exit nano:
# source ~/.bash_profile
### These are the locations in which we will put the important items below!!

# R  version
# mkdir -p ~/R-4.3.1/library
# echo 'R_LIBS_USER="~/R-4.3.1/library"' > $HOME/.Renviron
# export R_LIBS_USER="~/R-4.3.1"

## Get the github repo
# Follow these instructions on personal access tokens: https://chrisyeh96.github.io/2021/08/13/github-personal-access-tokens.html
# git clone https://github.com/flyconnectome/fafbpipeline.git
## Then log in with your github account
## To get a specific branch
# cd fafbpipeline
# git fetch origin flywire

# ssh-copy-id ab714@transfer.rc.hms.harvard.edu


# scp -r /Users/abates/.config/gcloud/legacy_credentials/alexander.shakeel.bates@gmail.com/adc.json ab714@transfer.rc.hms.harvard.edu:/home/ab714/.config/gcloud/legacy_credentials/alexander.shakeel.bates@gmail.com/
## Dropbox
# token <- drop_auth()
# saveRDS(token, "asbdropbox.RDS")
# scp -r /Users/abates/projects/flyconnectome/fafbpipeline/asbdropbox.RDS ab714@transfer.rc.hms.harvard.edu:/home/ab714/fafbpipeline/flywire/

## Install Java, from a local download:
## Get from :https://www.techspot.com/downloads/5553-java-jdk.html
# scp /Users/abates/Downloads/jdk-11.0.20_linux-x64_bin.tar.gz ab714@transfer.rc.hms.harvard.edu:/home/ab714/opt/jdk-11.0.20_linux-x64_bin.tar.gz
# ssh ab714@o2.hms.harvard.edu
# tar -xf opt/jdk-11.0.20_linux-x64_bin.tar.gz
# R CMD javareconf

## Install libgit2:
## Get from :https://github.com/libgit2/libgit2/releases
# scp /Users/abates/Downloads/libgit2-1.7.1.tar.gz ab714@transfer.rc.hms.harvard.edu:/home/ab714/opt/libgit2-1.7.1.tar.gz
# ssh ab714@o2.hms.harvard.edu
# tar -xf opt/libgit2-1.7.1.tar.gz

### blender, from a local download: https://www.blender.org/download/
# scp /Users/abates/Downloads/blender-3.6.1-linux-x64.tar.xz ab714@transfer.rc.hms.harvard.edu:/home/ab714/blender-3.6.1-linux-x64.tar.xz
# ssh ssh ab714@o2.hms.harvard.edu
# tar -xf blender-3.6.1-linux-x64.tar.xz

## Units
# wget ftp://ftp.unidata.ucar.edu/pub/udunits/udunits-2.2.28.tar.gz
# export INSTALL_PREFIX=$HOME/.local
# tar zxf udunits-2.2.28.tar.gz
# cd ./udunits-2.2.28/
# ./configure --enable-shared --prefix=$HOME/.local
# make
# make install
# export LD_LIBRARY_PATH=$INSTALL_PREFIX/lib:$LD_LIBRARY_PATH
# ls -l /home/ab714/.local/include/
# wget https://cran.r-project.org/src/contrib/Archive/udunits2/udunits2_0.13.tar.gz
# R CMD INSTALL udunits2_0.13.tar.gz --configure-args="--with-udunits2-include=$INSTALL_PREFIX/include --with-udunits2-lib=$INSTALL_PREFIX/lib/"
# wget https://cran.r-project.org/src/contrib/Archive/units/units_0.8-2.tar.gz
# R CMD INSTALL units_0.8-2.tar.gz --configure-args="--with-udunits2-include=$INSTALL_PREFIX/include --with-udunits2-lib=$INSTALL_PREFIX/lib/"

## download GDAL
# wget http://download.osgeo.org/gdal/2.2.3/gdal-2.2.3.tar.gz
# tar xzf gdal-2.2.3.tar.gz
# cd gdal-2.2.3
# ./configure --enable-shared --prefix=$HOME/.local
# make
# make install

## Install proj
# wget http://download.osgeo.org/proj/proj-7.0.0.tar.gz 
# tar zxvf proj-7.0.0.tar.gz 
# cd proj-7.0.0 \ 
#export PKG_CONFIG_PATH=/home/ab714/.local/lib:$PKG_CONFIG_PATH
# ./configure --enable-shared --prefix=/home/ab714/.local --datadir=/home/ab714/.local/share
#make 
#make install install-am 

## Google authentication, from a local copy:
# scp -r ~/Library/Caches/gargle/ ab714@transfer.rc.hms.harvard.edu:/home/ab714/Library/Caches/
# scp /Users/abates/Library/CloudStorage/GoogleDrive-ab2248@cam.ac.uk/Shared drives/hemibrainr/annotations/hemibrainr-571dd013f664.json ab714@transfer.rc.hms.harvard.edu:/fafbz/hemibrainr_data/annotations/

# scp -r /Users/abates/Library/CloudStorage/GoogleDrive-ab2248@cam.ac.uk/Shared\ drives/hemibrain/fafbsynapses/ ab714@transfer.rc.hms.harvard.edu:/n/data1/hms/neurobio/wilson/fafbz/fafbsynapses/
# scp -r /Users/abates/Library/CloudStorage/GoogleDrive-ab2248@cam.ac.uk/Shared\ drives/hemibrain/fafbsynapses/flywire_synapses.db ab714@transfer.rc.hms.harvard.edu:/n/data1/hms/neurobio/wilson/fafbz/fafbsynapses/flywire_synapses.db
# scp -r /Users/abates/Library/CloudStorage/GoogleDrive-ab2248@cam.ac.uk/Shared\ drives/hemibrain/fafbsynapses/synister_fafb_whole_volume_v3_t11.parquet ab714@transfer.rc.hms.harvard.edu:/n/data1/hms/neurobio/wilson/fafbz/fafbsynapses/synister_fafb_whole_volume_v3_t11.parquet
# scp -r /Users/abates/Library/CloudStorage/GoogleDrive-ab2248@cam.ac.uk/Shared\ drives/hemibrain/fafbsynapses/synister_fw_mat571_t11_synapses.feather ab714@transfer.rc.hms.harvard.edu:/n/data1/hms/neurobio/wilson/fafbz/fafbsynapses/synister_fw_mat571_t11_synapses.feather
# scp -r /Users/abates/Library/CloudStorage/GoogleDrive-ab2248@cam.ac.uk/Shared\ drives/hemibrain/fafbsynapses/synister_hemi_whole_volume_v0_t8.feather ab714@transfer.rc.hms.harvard.edu:/n/data1/hms/neurobio/wilson/fafbz/fafbsynapses/synister_hemi_whole_volume_v0_t8.feather

# rsync -r /Users/abates/Library/CloudStorage/GoogleDrive-ab2248@cam.ac.uk/Shared\ drives/hemibrain/fafbsynapses/synister_fafb_whole_volume_v3_t11.parquet ab714@transfer.rc.hms.harvard.edu:/n/data1/hms/neurobio/wilson/fafbz/fafbsynapses/synister_fafb_whole_volume_v3_t11.parquet
# rsync -r /home/ab714/rclone/flywire_neurons/flywire_data.sqlite /n/data1/hms/neurobio/wilson/fafbz/hemibrainr_data/flywire_neurons/flywire_data.sqlite
# rsync -r /home/ab714/rclone/hemibrain_neurons /n/data1/hms/neurobio/wilson/fafbz/hemibrainr_data/

# ## rclone
# curl -O https://downloads.rclone.org/rclone-current-linux-amd64.zip
# unzip rclone-current-linux-amd64.zip
# cd rclone-*-linux-amd64
# chmod 755 rclone
# rclone config
# --------------------
# Walk through `rclone config` interactively to authorise a new remote
# (the previous example config block has been redacted — it held OAuth
# tokens). The resulting `~/.config/rclone/rclone.conf` will look like:
#   [hemibrainr]
#   type = drive
#   client_id = <YOUR-CLIENT-ID>.apps.googleusercontent.com
#   client_secret = <YOUR-CLIENT-SECRET>
#   scope = drive
#   token = {...}                # populated by `rclone authorize`
#   team_drive = <YOUR-TEAM-DRIVE-ID>
#   root_folder_id =
# --------------------
# rclone mount hemibrainr rclone  

# # Install gettext
# wget https://ftp.gnu.org/gnu/gettext/gettext-0.22.tar.gz
# tar xvf gettext-0.22.tar.gz
# cd gettext-0.22
# ./configure --prefix=$HOME/.local
# make
# make install DESTDIR=$HOME/.local/usr/local/lib/
# export PKG_CONFIG_PATH=$HOME/.local/usr/local/lib/:$PKG_CONFIG_PATH
# export PATH=$HOME/.local/usr/local/lib/:$PATH

# # install fusermount3
# wget https://src.fedoraproject.org/repo/pkgs/fuse3/fuse-3.16.1.tar.gz/sha512/3f07919a7451a30d05bc174e2f8ec2c65b6225e63b4bbb40f2b097d760d4aa9b753a1da5da4874855094a01979fc4313ffabce54668ed20a6900f5eed92912c8/fuse-3.16.1.tar.gz
# signify -V -m fuse-3.16.1.tar.gz -p fuse-3.16.pub
# tar xvf fuse-3.16.1.tar.gz
# cd fuse-3.16.1
# meson setup ..
# DESTDIR=$HOME/.local ninja install
# ninja
# ninja install 

# Install Elastix
# wget https://github.com/SuperElastix/elastix/releases/download/5.1.0/elastix-5.1.0-linux.zip
# mkdir elastix-5.1.0-linux
# unzip elastix-5.1.0-linux.zip -d elastix-5.1.0-linux
# rm elastix-5.1.0-linux.zip

# gsutil cp gs://lee-lab_brain-and-nerve-cord-fly-connectome/m282/synapses_v1_id_size_prerootid_postrootid.parquet /n/data1/hms/neurobio/wilson/banc/
# gsutil cp gs://leelab_fly_cns/files/banc_nt_prediction_w_sizethresh_5_09072025.parquet /n/data1/hms/neurobio/wilson/banc/

## In R:
install.packages(c('usethis','remotes','googledrive','doMC', 'foreach','doParallel', 'reticulate', 'rJava', 'slackr', 'tidyverse','arrow', 'devtools', 'qs'))
remotes::install_github('natverse/nat.jrcbrains')
remotes::install_github('natverse/fafbseg')
remotes::install_github('flyconnectome/hemibrainr')
remotes::install_github('natverse/bancr')
library(nat.jrcbrains)
download_saalfeldlab_registrations()
library(fafbseg)
simple_python("full")
setwd("/n/data1/hms/neurobio/wilson/fafbz/")
hemibrainr:::hemibrainr_folder_structure("/n/data1/hms/neurobio/wilson/fafbz/hemibrainr_data")
dir.create("/n/data1/hms/neurobio/wilson/fafbz/hemibrainr_data/flywire_neurons/obj/")
dir.create("/n/data1/hms/neurobio/wilson/fafbz/hemibrainr_data/flywire_neurons/swc/")
dir.create("/n/data1/hms/neurobio/wilson/fafbz/hemibrainr_data/hemibrain_nblast/")
dir.create("/n/data1/hms/neurobio/wilson/fafbz/rtmp")

# Transfer in what data we already have
fw.neurons <- nat::neuronlist()
save(fw.neurons, file = "/n/data1/hms/neurobio/wilson/fafbz/hemibrainr_data/flywire_neurons/FlyWire/flywire_neurons.rda")

# rsync -r /home/ab714/rclone/flywire_neurons /n/data1/hms/neurobio/wilson/fafbz/hemibrainr_data/
# 
# 
# # Set up synapse table data bases
# srun --pty -p interactive -t 0-23:00 --mem 100G -c 1 bash
# module load gcc/9.2.0 R/4.3.1
# cd fafbpipeline/flywire
# R
# library(hemibrainr)
# library(DBI)
# fafbsynapses <- '/home/ab714/fafbsynapses/'
# fdb <- read.csv('/home/ab714/fafbsynapses/flywire_synapses.csv')
# con <- RSQLite::dbConnect(RSQLite::SQLite(), file.path(fafbsynapses,'flywire_synapses.db'), flags = RSQLite::SQLITE_RO)
# dbWriteTable('/home/ab714/fafbsynapses/flywire_synapses.db', "flywireids", fdb)

# 630 skeletons
cd /n/data1/hms/neurobio/wilson/fafbz/hemibrainr_data/flywire_neurons/FlyWire/630
wget https://zenodo.org/record/8077335/files/sk_lod1_630_healed_ds2.zip?download=1
unzip sk_lod1_630_healed_ds2.zip

# 783 skeletons
cd /n/data1/hms/neurobio/wilson/fafbz/hemibrainr_data/flywire_neurons/783
wget https://zenodo.org/records/10877326/files/sk_lod1_783_healed_ds2.parquet?download=1

# Install flywire environment
srun --pty -p interactive -t 0-10:00 --mem 5G -c 1 bash
module load miniconda3/4.10.3
module load python/3.13.1
conda create --name flywire_env python=3.13.1
source activate flywire_env
pip install --upgrade --prefer-binary numpy
pip install --upgrade --prefer-binary cloud-volume
pip install --upgrade --prefer-binary seatable_api!=2.6.3
pip install --upgrade --prefer-binary caveclient
pip install --upgrade --prefer-binary fafbseg
pip install --upgrade --prefer-binary skeletor
pip install --upgrade --prefer-binary fastremap
pip install --upgrade --prefer-binary ncollpyde
pip install --upgrade --prefer-binary meshparty
pip install --upgrade --prefer-binary pyembree

# Test
module load miniconda3/4.10.3
module load python/3.13.1
source activate flywire_env #/home/ab714/.local/share/r-miniconda/envs/r-reticulate
python
from cloudvolume import CloudVolume
import trimesh as tm
import skeletor as sk
m = tm.load_mesh('/n/data1/hms/neurobio/wilson/fafbz//hemibrainr_data/flywire_neurons//obj//720575940628524292.obj',process=False)
m = sk.utilities.make_trimesh(m, validate=False)
swc2 = sk.skeletonize.by_wavefront(mesh=m, waves=1, step_size=100, progress=False)
print(swc2)






# Get an O2 GPU
srun -n 1 --pty -t 1:00:00 -p gpu_quad --gres=gpu:rtx8000:1 --mem 5G bash
module load gcc/6.2.0 cuda/10.2 miniconda3/4.10.3 python/3.6.0

# Create a conda environment for deepCAD
conda init bash
conda create -n deepcad python=3.6
conda activate deepcad

# Install pytorch compatible with chosen CUDA
pip install torch==1.3.1
conda install pytorch==1.7.1 torchvision==0.8.2 torchaudio==0.7.2 cudatoolkit=10.2 -c pytorch

# Test cuda is available through torch
python
>>import torch
>>torch.cuda.is_available()
>>exit()

# Install other dependencies
conda install -c anaconda opencv scikit-learn scipy matplotlib
conda install scikit-image 
conda install -c conda-forge h5py pyyaml tensorboardx tifffile

# Clone deepCAD, this is a fork on the wilson-lab github
git clone https://github.com/wilson-lab/DeepCAD
  
# Download exemplar data
wget --no-check-certificate 'https://docs.google.com/uc?export=download&id=1k3de4r6fCs1kvpn0dcxSspJjZCVXbSnd' -O /home/ab714/DeepCAD/DeepCAD_pytorch/pth/ModelForPytorch/G_12_1199.pth
wget --no-check-certificate 'https://docs.google.com/uc?export=download&id=1E3uXGFQ8A8Z2wZmpMYLICuAhiMVpfdPr' -O /home/ab714/DeepCAD/DeepCAD_pytorch/pth/ModelForPytorch/para.yaml
wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=1WdRC-SaulA_CQfUWL3eSA1CEqBS9RjiV' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=1WdRC-SaulA_CQfUWL3eSA1CEqBS9RjiV" -O /home/ab714/DeepCAD/DeepCAD_pytorch/datasets/DataForPytorch/noisy_6000frames.tif && rm -rf /tmp/cookies.txt

# Example train
cd /home/ab714/DeepCAD/DeepCAD_pytorch/
conda activate deepcad
python script.py train

# Example test
conda activate deepcad
python script.py test

# Example debug mode, if you need it
python -m pdb test.py --denoise_model ModelForPytorch \
--datasets_folder DataForPytorch \
--test_datasize 6000

# If something with the environment goes wrong, wipe and start again:
conda remove -n deepcad --all
conda info --envs

# Run from scratch
screen -S deepcad
srun -n 1 --pty -t 26:00:00 -p gpu_quad --gres=gpu:rtx8000:1 --mem 24G bash
module load gcc/6.2.0 cuda/10.2 miniconda3/4.10.3 python/3.6.0
cd /home/ab714/DeepCAD/DeepCAD_pytorch/
conda activate deepcad
python script.py test


# Transfer your own data to the server




# DeepCAD_RT







# Get an O2 GPU
screen -S deepcadrt
srun -n 1 --pty -t 26:00:00 -p gpu_quad --gres=gpu:rtx8000:1,vram:24G  --mem 120G -c 4 bash
module load gcc/6.2.0 
module load cuda/10.2 
module load miniconda3/4.10.3 
module load python/3.6.0

# # Install cudNN
# rsync -r /Users/abates/Downloads/cudnn-linux-x86_64-8.7.0.84_cuda10-archive.tar.xz ab714@transfer.rc.hms.harvard.edu:/home/ab714
# tar -xf cudnn-linux-x86_64-8.7.0.84_cuda10-archive.tar.xz
# cp cudnn-linux-x86_64-8.7.0.84_cuda10-archive/include/cudnn*.h /home/ab714/.local/include 
# cp -P cudnn-linux-x86_64-8.7.0.84_cuda10-archive/lib/libcudnn* /home/ab714/.local/lib 
# chmod a+r /home/ab714/.local/include/cudnn*.h /home/ab714/.local/lib/libcudnn*
#   
# # Add cuDNN to path
# export CPATH=$CPATH:$HOME/.local/include 
# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/.local/lib

# Create a conda environment for deepCAD
conda init bash
conda create -n deepcadrt python=3.6
source activate deepcadrt

# Install pytorch compatible with chosen CUDA
pip install torch==1.8.0
conda install pytorch==1.8.0 torchvision==0.8.2 torchaudio==0.7.2 cudatoolkit=10.2 -c pytorch
pip install --upgrade pip setuptools wheel
pip3 install opencv-python==4.1.2.30
pip install deepcad

## Install repo
git clone https://github.com/wilson-lab/DeepCAD-RT

# Try out
conda activate deepcadrt
cd /home/ab714/DeepCAD-RT/DeepCAD_RT_pytorch/
python demo_test_pipeline.py



# Movew data off of server to /files
screen -S deepcad
rsync -hvrPt /n/data1/hms/neurobio/wilson/DeepCAD_datasets/ /files/Neurobio/wilsonlab/DeepCAD_datasets/



  
#   file <- "/Users/abates/Library/CloudStorage/GoogleDrive-ab2248@cam.ac.uk/Shared drives/hemibrain/fafbsynapses/flywire_synapses.db"
# con <- RSQLite::dbConnect(RSQLite::SQLite(), db=file, flags = RSQLite::SQLITE_RO)
# RSQLite::dbListTables(con)
# res <- dplyr::tbl(con, 'flywireids')
# csv <- res %>%
#   dplyr::collect() %>%
#   as.data.frame()
# readr::write_csv(x=csv,
#                  file="/Users/abates/Library/CloudStorage/GoogleDrive-ab2248@cam.ac.uk/Shared drives/hemibrain/fafbsynapses/flywire_synapses.csv")
# 
# file <- "/Users/abates/Library/CloudStorage/GoogleDrive-ab2248@cam.ac.uk/Shared drives/hemibrain/fafbsynapses/20191211_fafbv14_buhmann2019_li20190805_nt20201223.db"
# con <- RSQLite::dbConnect(RSQLite::SQLite(), db=file, flags = RSQLite::SQLITE_RO)
# RSQLite::dbListTables(con)
# res <- dplyr::tbl(con, 'synlinks')
# csv <- res %>%
#   dplyr::collect() %>%
#   as.data.frame()
# readr::write_csv(x=csv,
#                  file="/Users/abates/Library/CloudStorage/GoogleDrive-ab2248@cam.ac.uk/Shared drives/hemibrain/fafbsynapses/20191211_fafbv14_buhmann2019_li20190805_nt20201223.csv")
# 
# 
# rsync -r /Users/abates/Library/CloudStorage/GoogleDrive-ab2248@cam.ac.uk/Shared\ drives/hemibrain/fafbsynapses/flywire_synapses.csv ab714@transfer.rc.hms.harvard.edu:/home/ab714/fafbsynapses/flywire_synapses.csv
# rsync -r /Users/abates/Library/CloudStorage/GoogleDrive-ab2248@cam.ac.uk/Shared\ drives/hemibrain/fafbsynapses/20191211_fafbv14_buhmann2019_li20190805_nt20201223.csv ab714@transfer.rc.hms.harvard.edu:/home/ab714/fafbsynapses/f20191211_fafbv14_buhmann2019_li20190805_nt20201223.csv
# 
# nohup scp -r /Users/abates/Library/CloudStorage/GoogleDrive-ab2248@cam.ac.uk/Shared\ drives/hemibrain/fafbsynapses/flywire_synapses.csv ab714@transfer.rc.hms.harvard.edu:/home/ab714/fafbsynapses/flywire_synapses.csv
# nohup rsync -r /Users/abates/Library/CloudStorage/GoogleDrive-ab2248@cam.ac.uk/Shared\ drives/hemibrain/fafbsynapses/20191211_fafbv14_buhmann2019_li20190805_nt20201223.csv ab714@transfer.rc.hms.harvard.edu:/home/ab714/fafbsynapses/f20191211_fafbv14_buhmann2019_li20190805_nt20201223.csv



# JOBID PARTITION     USER  ACCOUNT   PRIORITY       SITE        AGE      ASSOC  FAIRSHARE    JOBSIZE  PARTITION    QOSNAME        QOS        NICE                 TRES
# 16578366 gpu_quad     ab714   wilson     644382          0       5973          0     295540         13     342857 gpuquad_qo          0           0                     

# Account                    User  RawShares  NormShares    RawUsage  EffectvUsage  FairShare 
# -------------------- ---------- ---------- ----------- ----------- ------------- ---------- 
#   wilson                    ab714          1    0.000277     1227203      0.002920   0.295540 





srun --pty -p interactive -t 0-2:00 --mem 6G -c 1 bash
gsutil cp /n/data1/hms/neurobio/wilson/banc/connectivity/banc_data.sqlite gs://lee-lab_brain-and-nerve-cord-fly-connectome/v282
gsutil cp /n/data1/hms/neurobio/wilson/banc/meta/banc_manc_v1.2.1_nblast.csv gs://lee-lab_brain-and-nerve-cord-fly-connectome/nblast
gsutil cp /n/data1/hms/neurobio/wilson/banc/meta/banc_hemibrain_v1.2.1_nblast.csv gs://lee-lab_brain-and-nerve-cord-fly-connectome/nblast
gsutil cp /n/data1/hms/neurobio/wilson/banc/meta/banc_fafb_783_nblast.csv gs://lee-lab_brain-and-nerve-cord-fly-connectome/nblast
gsutil cp /n/data1/hms/neurobio/wilson/banc/meta/banc_mirror_nblast.csv gs://lee-lab_brain-and-nerve-cord-fly-connectome/nblast

gsutil cp /n/data1/hms/neurobio/wilson/banc/meta/2024-09-02_community_challenge_abdominal_neuromere.csv gs://lee-lab_brain-and-nerve-cord-fly-connectome/neuron_searches
gsutil cp /n/data1/hms/neurobio/wilson/banc/meta/2024-09-02_manc_efferent_matches.csv gs://lee-lab_brain-and-nerve-cord-fly-connectome/neuron_searches
gsutil cp /n/data1/hms/neurobio/wilson/banc/meta/2024-09-02_manc_sensory_matches.csv gs://lee-lab_brain-and-nerve-cord-fly-connectome/neuron_searches

gsutil ls gs://lee-lab_brain-and-nerve-cord-fly-connectome/imported_meshes/manc_v1.2.1_meshes_elastix_tpsreg_240721/meshes/** | wc -l
gsutil ls gs://lee-lab_brain-and-nerve-cord-fly-connectome/imported_meshes/fafb_783_meshes_elastix_tpsreg_240721/meshes/** | wc -l
gsutil ls gs://lee-lab_brain-and-nerve-cord-fly-connectome/imported_meshes/manc_v1.2.1_meshes_elastix_tpsreg_240721/meshes/** | wc -l

gsutil cp manc_banc_space_swc.zip gs://lee-lab_brain-and-nerve-cord-fly-connectome/skeletons

gsutil ls gs://dkronauer-ant-001-segmentations-prod/nucleus/240910-nucleus/** | wc -l


mv /n/data1/hms/neurobio/wilson/banc/matching/hemibrain/elastix_tpsreg_240721/CX_columnarby_cell_type/ /n/data1/hms/neurobio/wilson/banc/matching/hemibrain/elastix_tpsreg_240721/CX_columnar_by_cell_type/