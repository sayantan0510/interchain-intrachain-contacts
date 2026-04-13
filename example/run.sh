# for intra_chain_contact
python ../scripts/intra_chain_contact.py -s ./input/system.tpr -f ./input/system.xtc --out-prefix ./results/intra-chain_bb --contact-mode bb --chain-by auto --r0-nm 0.6 --tmin-ps 0.0 --tmax-ps 400000.0 --stride 1 --verbose --plot-chainwise-hists --write-chainwise-csv

# for inter_chain_contact
python ../scripts/inter_chain_contact.py -s ./input/system.tpr -f ./input/system.xtc --out-prefix ./results/INTER-chain_bb --contact-mode bb --chain-by auto --r0-nm 0.6 --tmin-ps 0.0 --tmax-ps 400000.0 --stride 1 --verbose --plot-chainwise-hists --write-chainwise-csv
