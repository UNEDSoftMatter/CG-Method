# /*
#  * Filename   : preprocess.sh
#  *
#  * Created    : 26.04.2016
#  *
#  * Modified   : mié 27 abr 2016 18:40:01 CEST
#  *
#  * Author     : jatorre
#  *
#  * Run        :  
#  *
#  * Purpose    :  
#  *
#  */


# let NSTEPS=2000
# 
# for ((i=1;i<=NSTEPS;i++));
# do
#   split --number=$i/NSTEPS output.positions | tail -n +10 | sort -n |awk '{print $2,$3,$4,$5}' > ./data/$i.pos
#   split --number=$i/NSTEPS output.velocities | tail -n +10 | sort -n |awk '{print $3,$4,$5}'    > ./data/$i.vel
#   ./CG $i > ./log/$i.log
#   rm ./data/$i.pos ./data/$i.vel
# done

for i in $(ls ./data/positions |grep .pos); do
  var=$(basename $i .pos)
  ./CG ${var} > ./log/${var}.log 
done
