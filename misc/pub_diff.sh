cat ../database/vdjdb-2018-01-17/vdjdb.txt | cut -f13 | sort | uniq > pub_prev
cat ../database/vdjdb-2018-05-05/vdjdb.txt | cut -f13 | sort | uniq > pub_new
diff -a pub_prev pub_new