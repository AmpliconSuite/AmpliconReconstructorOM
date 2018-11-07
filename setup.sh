#!/usr/bin/env bash

#This sets up AmpliconReconstuctor and SegAligner as a local command line utility, callable from any directory

#add AmpliconReconstructor to path
echo "export AR_SRC=$PWD" >> ~/.bashrc
echo "alias AmpliconReconstructor='$AR_SRC/AmpliconReconstructor'" >> ~/.bashrc
echo "alias AR='$AR_SRC/AmpliconReconstructor'" >> ~/.bashrc

#make SegAligner
cd SegAligner
make clean; make

#add SegAligner to path
echo "export SA_SRC=$PWD" >> ~/.bashrc
echo "alias SegAligner='$SA_SRC/SegAligner'"
