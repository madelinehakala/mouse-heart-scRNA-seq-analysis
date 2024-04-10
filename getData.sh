mkdir -p dataZone2
mkdir -p dataZone3

wget -O dataZone2/matrix.mtx.gz https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3885nnn/GSM3885059/suppl/GSM3885059_AVNmatrix.mtx.gz
wget -O dataZone2/barcodes.tsv.gz https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3885nnn/GSM3885059/suppl/GSM3885059_AVNbarcodes.tsv.gz
wget -O dataZone2/features.tsv.gz https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3885nnn/GSM3885059/suppl/GSM3885059_AVNgenes.tsv.gz

wget -O dataZone3/matrix.mtx.gz https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3885nnn/GSM3885060/suppl/GSM3885060_LPFmatrix.mtx.gz
wget -O dataZone3/barcodes.tsv.gz https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3885nnn/GSM3885060/suppl/GSM3885060_LPFbarcodes.tsv.gz
wget -O dataZone3/features.tsv.gz https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3885nnn/GSM3885060/suppl/GSM3885060_LPFgenes.tsv.gz
