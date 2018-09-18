
# python3
# Adrian Garcia Moreno
# for file in $(ls -1); do awk '{print "chr"$0}' $file >"chr"$file ; done

indir = "rdm_hyi/"
out_dir = "RandomModels"
cuts = " ".join(map(str, [4.5,5]))

import os

rdmmodels = [indir+rdmmodel for rdmmodel in os.listdir(indir)]
for rdmmodel in rdmmodels:
    cmd = "python3 criteria_divider.py -file {} -f_clm 4 -s_clm 5 -cuts {} -o {} &".format(rdmmodel, cuts, out_dir)
    print(cmd)
    os.system(cmd)

toannotedirs = ["{}/{}/Cdivided/".format(out_dir,toannotedir) for toannotedir in os.listdir(out_dir)]
toenrdirs = []
for toannotedir in toannotedirs:
    basedir = os.path.dirname(os.path.dirname(toannotedir))
    outdir = "{}/GOenr".format(basedir)
    toenrdirs.append(outdir+"/")
    cmd = "python3 bedtooling.py -d {} -annot annotations/GeneIDs_interim_GRCh37.p13_top_level_2017-01-13.gff3.gz -o {}".format(toannotedir, outdir)

for toenrdir in toenrdirs:
    basedir = os.path.dirname(os.path.dirname(toenrdirs[0]))
    outdir = "{}/GOenr".format(basedir)
    cmd = "Rscript TopGOer.r {} annotations/9606_geneID2GO.map 10 {} -mode CC,BP -pval_thres 0.001 ".format(indir, outdir)
