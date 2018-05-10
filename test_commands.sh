### criteria_divider:

## File 
# Criteria + CutOff
python3 criteria_divider.py -file example_data/chrs+feature+score/1.bed -f_clm 4 -s_clm 5 -cut 2 -o example_data/example_results/FCC/
# Criteria
python3 criteria_divider.py -file example_data/chrs+feature+score/1.bed -f_clm 4 -o example_data/example_results/FC/
# CutOff
python3 criteria_divider.py -file example_data/chrs+feature+score/1.bed -s_clm 5 -cut 3 -o example_data/example_results/FC-O/

## Directory
# Criteria + CutOff
python3 criteria_divider.py -dir example_data/chrs+feature+score/ -f_clm 4 -s_clm 5 -cut 3 -o example_data/example_results/DCC/
# Criteria
python3 criteria_divider.py -dir example_data/chrs+feature+score/ -f_clm 4 -o example_data/example_results/DC/
# CutOff
python3 criteria_divider.py -dir example_data/chrs+feature+score/ -s_clm 5 -cut 3 -o example_data/example_results/DC-O/


### gff3_manager

## Download Without Transform
python3 gff3_manager.py -down_gff3 gcode_hg19
python3 gff3_manager.py -down_gff3 gcode_hg20
python3 gff3_manager.py -down_gff3 ncbi_hg19
python3 gff3_manager.py -down_gff3 ncbi_hg20

## Download With Transform
python3 gff3_manager.py -down_gff3 ncbi_hg19 -to_chr
python3 gff3_manager.py -down_gff3 ncbi_hg20 -to_chr

## File as input and Transform
python3 gff3_manager.py -gff3 interim_GRCh37.p13_top_level_2017-01-13.gff3.gz -consortium ncbi -to_chr
python3 gff3_manager.py -gff3 ref_GRCh38.p12_top_level.gff3.gz -consortium ncbi -to_chr

## File without transform
python3 gff3_manager.py -gff3 gencode.v28.annotation.gff3.gz -consortium gcode
