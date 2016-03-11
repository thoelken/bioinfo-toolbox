#!/bin/bash

gawk '!/circRNA_start/ {print $2"\t"$3"\t"$4"\tciri:"$5":"$7"\t0\t"$11}' $1
