#! /usr/bin/bash

gcc -c -g read_config_para.c
gcc -c -g read_interface_material.c
gcc -c -g pre_para_assign.c
gcc -c -g pre_cal_dip_area.c
gcc -c -g pre_media_parameterization.c 
gcc -c -g FD2D_ModelPrepare.c

gcc -o pre FD2D_ModelPrepare.o read_config_para.o read_interface_material.c pre_para_assign.o pre_cal_dip_area.o pre_media_parameterization.o -lm
#gcc -o b main.o read_conf.o

rm *.o
