awk -F\t ' {if(match($2,", ")) {split($2,a,", ") ; for (i=1; i<=length(a); i++) {sub($2, a[i]) ; print}} else print}'
