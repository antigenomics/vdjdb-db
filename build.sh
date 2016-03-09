mkdir database/
cd src/
python build_db.py
Rscript db_summary.R | tee ../database/vdjdb.summary.txt
cd ..

STAMP=`date +%Y-%m-%d`-`git rev-parse --short HEAD`

cd database/
cp ../metadata/vdjdb.meta.txt vdjdb.meta.txt
cp ../LICENSE.txt LICENSE.txt
zip -r vdjdb-$STAMP.zip vdjdb.meta.txt vdjdb.txt LICENSE.txt