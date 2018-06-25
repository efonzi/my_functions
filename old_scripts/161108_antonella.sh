#!/bin/bash

dos2unix -c Mac |
awk -F\t 'NR>1 {OFS="\t"}{sub(/"/, "", $2)};{print}' |
awk -F\t 'NR>1 {OFS="\t"}{sub(/"/, "", $2)};{print}' |
awk -F\t '{if(match($2,", ")) {split($2,a,", ") ; for (i=1; i<=length(a); i++) {sub($2, a[i]) ; print}} else print}' |
awk -F\t 'NR==1 {print}; NR>1 {print | "sort -u"}'

#awk -F\t '{NR==1 print $0}{print $0 | "sort -u"}'
