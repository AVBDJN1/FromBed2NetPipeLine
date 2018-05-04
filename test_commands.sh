### CriteriaDivider:

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
