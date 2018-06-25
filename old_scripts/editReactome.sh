dos2unix -c Mac |
awk -F\t '{OFS="\t"}{sub(/NA/, ".", $2)};{print}' |
awk -F\t '{OFS="\t"}{sub(/PATH.ME/, "PATHNAME", $2)};{print}' |
awk -F\t '{OFS="\t"}{sub(/.DE/, "NADE", $2)};{print}' |
awk -F\t 'NR>1 {OFS="\t"}{sub(/NA/, ".", $3)};{print}' |
awk -F\t 'NR>1 {OFS="\t"}{sub(/Homo sapiens: /, "", $0)};{print}'
