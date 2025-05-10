# Compile database

cd py_src/
python runBuidDatabase.py
cd ..

# Summary HTML

cd summary/
Rscript -e "rmarkdown::render('vdjdb_summary.Rmd')"
python MakeEmbedableHtml.py
cd ..


# Compute motifs
if [ -d "../vdjdb-motifs" ] 
then
  echo "Found 'vdjdb-motifs' folder. Processing..."
  
  mkdir ../vdjdb-motifs/vdjdb_dump/
  cp -r database/*.txt ../vdjdb-motifs/vdjdb_dump/

  cd ../vdjdb-motifs/
  Rscript -e "rmarkdown::render('compute_vdjdb_motifs.Rmd')"
  cp cluster_members.txt ../vdjdb-db/database/
  cp motif_pwms.txt ../vdjdb-db/database/
  cd ../vdjdb-db/
else
  echo "'vdjdb-motifs' folder is missing. Skipping..."
fi

cd py_src/
python MotifsScoresAssembler.py
cd ..

# Gather database_

ls ../vdjdb-motifs/
echo '_____'
ls ../vdjdb-motifs/vdjdb_export
echo '_____'
ls ../vdjdb-motifs/tcrnet


cd database/
DD=`date +%F`
mkdir vdjdb-$DD
cp ../summary/vdjdb_summary_embed.html vdjdb-$DD/
cp *.txt vdjdb-$DD/

# Update latest version

RR="https://github.com/antigenomics/vdjdb-db/releases/download/$DD/vdjdb-$DD.zip"
RRP=`cat ../latest-version.txt | head -n 1`

if [[ $RRP != $RR ]];
then
  echo "Releasing"
  echo $RR  | cat - ../latest-version.txt > temp && mv temp ../latest-version.txt
else
  echo "Overwriting"
fi

cp ../LICENSE.txt vdjdb-$DD/
cp ../latest-version.txt vdjdb-$DD/

# ZIP together

zip -rj vdjdb-$DD.zip vdjdb-$DD/

echo "scp `hostname`:`pwd`/vdjdb-$DD.zip ."