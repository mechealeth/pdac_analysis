import os, sys
import pandas as pd
from pathlib import Path

os.chdir("/lustre/home/whhou/00.datasets/2021_Data/process/mapping/PDAC/ST")

# raw_meta = pd.read_csv("/lustre/home/whhou/00.datasets/2021_Data/process/data/PDAC/ST/PDAC_ST_20220522.csv")
# raw_meta = pd.read_csv("/lustre/home/whhou/00.datasets/2021_Data/process/data/PDAC/ST/PDAC_FF_202208.csv")
raw_meta = pd.read_csv("/lustre/home/whhou/00.datasets/2021_Data/process/data/PDAC/ST/PDAC_FF_202209.csv")


#only keep one record
raw_meta.drop_duplicates(subset="SampleName", keep="first", inplace=True)
raw_meta = raw_meta.reset_index()

for l in range(0, raw_meta.shape[0]):
	if not os.path.exists(raw_meta.SampleName[l]):
		os.makedirs(raw_meta.SampleName[l])

	# lurmTmp = "/lustre/home/whhou/00.datasets/2021_Data/process/mapping/PDAC/ST/spaceranger_count_FFPE_tmp.slurm"
	lurmTmp = "/lustre/home/whhou/00.datasets/2021_Data/process/mapping/PDAC/ST/spaceranger_count_FF_tmp.slurm"

	sampleslurm = os.path.join(raw_meta.SampleName[l], raw_meta.SampleName[l]+'.slurm')
	os.system('cp %s %s'%(lurmTmp, sampleslurm))

	FQdir = os.path.dirname(raw_meta.targetFile[l])
	FQdir = FQdir.replace('/', '\/')
	cmd = "sed -i 's/FQDIR/%s/g' %s"%(FQdir, sampleslurm)
	os.system(cmd)

	cmd = "sed -i 's/SAMPLENAME/%s/g' %s"%(raw_meta.SampleName[l], sampleslurm)
	os.system(cmd)

	imgFILE = raw_meta.imgFile[l]
	imgFILE = imgFILE.replace('/', '\/')
	cmd = "sed -i 's/IMAGEFILE/%s/g' %s"%(imgFILE, sampleslurm)
	os.system(cmd)

	SLIDEFILE = raw_meta.Slide[l]
	SLIDEFILE = SLIDEFILE.replace('/', '\/')
	cmd = "sed -i 's/SLIDEFILE/%s/g' %s"%(SLIDEFILE, sampleslurm)
	os.system(cmd)

	SLIDE = Path(raw_meta.Slide[l]).stem
	cmd = "sed -i 's/SLIDE/%s/g' %s"%(SLIDE, sampleslurm)
	os.system(cmd)


	AREAINDEX = raw_meta.areaIndex[l]
	cmd = "sed -i 's/AREA/%s/g' %s"%(AREAINDEX, sampleslurm)
	os.system(cmd)

	print('successfully prepared %s'%sampleslurm)
