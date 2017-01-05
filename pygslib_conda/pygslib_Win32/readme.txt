From http://stackoverflow.com/questions/33709391/using-multiple-python-engines-32bit-64bit-and-2-7-3-5

set CONDA_FORCE_32BIT=1
conda create -n py27_32 python=2.7

set CONDA_FORCE_32BIT=1
activate py27_32

#Deactivate it:
deactivate py27_32