#!/usr/bin/env bash

#This sets up AmpliconReconstuctor and SegAligner as a local command line utility, callable from any directory

#add AmpliconReconstructor and SegAligner to path
echo "export AR_SRC=$PWD" >> ~/.bashrc
echo "export SA_SRC=$PWD" >> ~/.bashrc
source ~/.bashrc

#alias AmpliconReconstructor
echo "alias AmpliconReconstructor='python2 $AR_SRC/AmpliconReconstructor.py'" >> ~/.bashrc
echo "alias AR='python2 $AR_SRC/AmpliconReconstructor.py'" >> ~/.bashrc

#make SegAligner
cd SegAligner
make clean; make

#alias SegAligner
echo "alias SegAligner='$SA_SRC/SegAligner'" >> ~/.bashrc

#make changes effective for current session
source ~/.bashrc

echo -e "\nSetup complete\n"
