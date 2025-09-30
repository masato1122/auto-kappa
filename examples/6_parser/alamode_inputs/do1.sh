dir1=../../../../../ALL/2500-5000/mp-3614
#dir1=../../../../../ALL/2500-5000/mp-4324

cp $dir1/relax/CONTCAR ./FILES/POSCAR-primitive
cp $dir1/relax/freeze-1/CONTCAR ./FILES/POSCAR-primitive
cp $dir1/harm/force/prist/POSCAR ./FILES/POSCAR-supercell
cp $dir1/result/DFSET.cube_fd ./FILES/DFSET.cube
cp $dir1/result/FC2.xml ./FILES/
cp $dir1/result/FC3_fd.xml ./FILES/
cp $dir1/harm/bandos/BORNINFO ./FILES/

