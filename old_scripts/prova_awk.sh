awk 'NR>7' |
awk -F\t '{OFS="\t"}{for(i=1; i<=NF; i++) if($i ~ /^ *$/) $i = "."};{print}' |
awk -F\t '{OFS="\t"}{if(match($15,"mi")) {} else $15="."; print}' |
awk '{gsub(/ - /, "_")};{print}' |
awk '{sub(/Region/, "chr\tstart\tend")};{print}' |
awk '{gsub(/Gene Symbols/, "Gene Symbol")};{print}' |
awk -F\t 'NR>1 {OFS="\t"}{sub(/:/, "\t", $1)};{print}' |
awk -F\t 'NR>1 {OFS="\t"}{sub(/,/, "", $2)};{print}' |
awk -F\t 'NR>1 {OFS="\t"}{sub(/,/, "", $2)};{print}' |
awk -F\t 'NR>1 {OFS="\t"}{sub(/-/, "\t", $2)};{print}' |
awk -F\t 'NR>1 {OFS="\t"}{sub(/,/, "", $3)};{print}' |
awk -F\t 'NR>1 {OFS="\t"}{sub(/,/, "", $3)};{print}' |
awk -F\t ' {if(match($15,", ")) {split($15,a,", ") ; for (i=1; i<=length(a); i++) {sub($15, a[i]) ; print}} else print}'
