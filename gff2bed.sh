#!/bin/bash

perl -ne 's/^(\S+)\t\S+\t\S+\t(\d+)\t(\d+)\t\S+\t([-\+])\t\S+\tID=([^;]+);?(.*)$/\1\t\2\t\3\t\5\t0\t\4/ && print' $1
