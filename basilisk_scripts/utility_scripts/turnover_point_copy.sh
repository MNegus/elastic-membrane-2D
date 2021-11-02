

for DIRNAME in dirichlet_surface neumann_surface no_surface
do
    for MAXLEVEL in 10 11 12 13
    do
        ORIG=/scratch/negus/jet_root_height_validation/$DIRNAME/max_level_$MAXLEVEL/raw_data/
        DEST=~/Desktop/jet_energy/$DIRNAME/max_level_$MAXLEVEL/
        cp $ORIG/turnover_points_basilisk.txt $DEST
    done
done
