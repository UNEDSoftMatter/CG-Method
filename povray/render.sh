for i in $(ls |grep x); do

  sed 's/XXXXXX/'$i'/' temperature.pov > $i.pov

povray +I$i.pov +O$i.png \
       +W800 +H450 \
       +FN16 +Q11  \
       +A0.01 +AM2 \
       +J0.5       \
       -D

done

#       +W1920 +H1080 \

# povray +I$1.pov +O$1.png \
#        +W800 +H600 \
#        +FN16 +Q11 \
#        +A0.01 +AM2 \
#        +J0.5 
# 
