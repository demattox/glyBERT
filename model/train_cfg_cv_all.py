
import subprocess


lectinlist=['GSL','H3N8','UEAI','LCA','SNA','PNA','RCAI','SBA','PHA-E','PHA-L','WGA','jacalin','ConA','MAL_II','MAL-I','ABA','HA','DBA','PSA','Human']


for i in range(1,6):
    for l in lectinlist:
        subprocess.call("python train_tune_cv.py "+l+'_'+str(i), shell=True)