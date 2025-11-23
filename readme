for i in insar/???A* insar/???D* ; do echo $i >> frames.txt ; done
cat frames.txt | awk '{print "insardir_"NR":   "$0"/"}'
rm -rf frames.txt
