Rscript preprocess.R
mvn clean install
java -Xmx4G -jar target/cdr3net-0.0.1.jar 3 0 vdjdb.nodes.txt vdjdb.edges.txt
