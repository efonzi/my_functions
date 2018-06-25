awk '/Missense/{gsub($7, "missense/nonsense")};{print}' |
awk '/Unknown/{gsub($7, "unknown")};{print}' |
awk '/Synonymous/{gsub($7, "silent")};{print}' |
awk '/Nonsense/{gsub($7, "missense/nonsense")};{print}' |
awk '/Insertion/{gsub($7, "insertion/deletion")};{print}' |
awk '/Deletion/{gsub($7, "insertion/deletion")};{print}' |
awk '/Frameshift/{gsub($7, "insertion/deletion")};{print}' |
awk '/Silent/{gsub($7, "silent")};{print}' |
awk '/Noncoding/{gsub($7, "noncoding")};{print}' |
awk '/Splicing/{gsub($7, "splicing")};{print}'
