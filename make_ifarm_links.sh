#!/usr/bin/bash
# This script generates useful links (to DATA, ROOT and other directories) for ifarm machines
ln -s /volatile/halld/home/nseptian/TRDROOT ROOT
ln -s /volatile/halld/home/nseptian/TRDRootOutput TRDRootOutput
mkdir RootOutput


# ln -s /home/hdtrdops/DATA DATA
# ln -s /gluonraid3/data4/rawdata/trd/LOG LOG

# ln -s /home/hdtrdops/soft_ml4fpga/setup_env.csh setup_env.csh
# ln -s /home/hdtrdops/soft_ml4fpga/setup_env.sh setup_env.sh
# mkdir mlpOutput