#!/bin/bash

dos2unix -c Mac |
awk '{sub(/Chromosome Region/, "chr\tstart\tend")};{print}' |
awk -F\t 'NR>1 {OFS="\t"}{sub(/"/, "", $2)};{print}' |
awk -F\t 'NR>1 {OFS="\t"}{sub(/:/, "\t", $2)};{print}' |
awk -F\t 'NR>1 {OFS="\t"}{sub(/,/, "", $3)};{print}' |
awk -F\t 'NR>1 {OFS="\t"}{sub(/,/, "", $3)};{print}' |
awk -F\t 'NR>1 {OFS="\t"}{sub(/-/, "\t", $3)};{print}' |
awk -F\t 'NR>1 {OFS="\t"}{sub(/,/, "", $4)};{print}' |
awk -F\t 'NR>1 {OFS="\t"}{sub(/,/, "", $4)};{print}' |
awk -F\t 'NR>1 {OFS="\t"}{sub(/"/, "", $4)};{print}' |
awk -F\t 'NR>1 {OFS="\t"}{sub(/"/, "", $12)};{print}' |
awk -F\t 'NR>1 {OFS="\t"}{sub(/"/, "", $12)};{print}' |
awk -F\t '{if(match($12,", ")) {split($12,a,", ") ; for (i=1; i<=length(a); i++) {sub($12, a[i]) ; print}} else print}'
