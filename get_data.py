import os
import gzip

wget_GEO = "wget 'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3885nnn/GSM3885059/suppl/GSM3885059%5FAVNmatrix.mtx.gz'"
os.system(wget_GEO)

with gzip.open('GSM3885059_AVNmatrix.mtx.gz', 'rb') as f_in:
  with open('ZoneII.mtx', 'wb') as f_out:
    f_out.write(f_in.read())

wget_GEO = "wget 'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3885nnn/GSM3885060/suppl/GSM3885060%5FLPFmatrix.mtx.gz'"
os.system(wget_GEO)

with gzip.open('GSM3885059_AVNmatrix.mtx.gz', 'rb') as f_in:
  with open('ZoneIII.mtx', 'wb') as f_out:
    f_out.write(f_in.read())
