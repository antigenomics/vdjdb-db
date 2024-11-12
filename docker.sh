# /bin/sh'
mkdir -p vcs
cd vcs
echo $(pwd)
git clone https://github.com/antigenomics/vdjdb-db vdjdb-db
git clone https://github.com/antigenomics/vdjdb-motifs vdjdb-motifs
git clone https://github.com/antigenomics/mirpy.git
cd vdjdb-db
git checkout pyVDJdb
echo $(pwd)
mkdir -p /root/output
bash release.sh 2>&1 | tee /root/output/buildlog
cp -r database/*zip /root/output/