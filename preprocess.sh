# /*
#  * Filename   : preprocess.sh
#  *
#  * Created    : 26.04.2016
#  *
#  * Modified   : mar 26 abr 2016 16:13:52 CEST
#  *
#  * Author     : jatorre
#  *
#  * Run        :  
#  *
#  * Purpose    :  
#  *
#  */


let NSTEPS=100

for ((i=1;i<=NSTEPS;i++));
do
  split --number=$i/NSTEPS $1 | tail -n +10 | sort -n |awk '{print $2,$3,$4,$5}' > ./data/$i.pos
  split --number=$i/NSTEPS $1 | tail -n +10 | sort -n |awk '{print $3,$4,$5}'    > ./data/$i.vel
  ./CG $i > ./log/$i.log
  rm ./data/$i.pos ./data/$i.vel
done
