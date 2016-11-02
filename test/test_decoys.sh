#!/bin/bash

# Top-level directory of package
PACKAGE_DIR=".."
PDBID="4ah2"
R_CH="A"
L_CH="C"
N_WINDOWS=3

complexid=$PDBID$R_CH$L_CH

# Generate pairwise decoys
python generate_decoys.py
# Create database
python $PACKAGE_DIR/scripts/create_database.py -w . -p $PDBID -d "${PDBID}_data.csv"
mv "scores_$PDBID.db" "scores_$complexid.db"
# Load decoys into database
python $PACKAGE_DIR/scripts/load_model_scores.py -d "scores_$complexid.db" -p $PDBID -r $R_CH -l $L_CH
# Perform combinatorial search
python $PACKAGE_DIR/scripts/find_paths_stepwise.py $complexid -d . -n $N_WINDOWS
# Compute occupancy score
python $PACKAGE_DIR/scripts/compute_occupancy_score.py $PDBID -r $R_CH -l $L_CH -n $N_WINDOWS -d .
# Rank and select paths
python $PACKAGE_DIR/scripts/select_paths.py $PDBID -r $R_CH -l $L_CH -n $N_WINDOWS -d .
