#!/usr/bin/bash

gcc -c -g Elastic2d_tti.c
gcc -c -g pre_read_conf.c

gcc -o b *.o
