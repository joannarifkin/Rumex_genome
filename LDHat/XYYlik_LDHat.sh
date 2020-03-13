for LG in 7 10 3 5 8
do

rm new_lk.txt rates.txt bounds.txt res.txt 
/ohta/felix.beaudry/scripts/LDhat-master/complete -n 54 -theta 0.006479939 -rhomax 100 -n_pts 101 
cp new_lk.txt lk_XYY.LG${LG}.LDhat.txt

done