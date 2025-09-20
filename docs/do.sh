
make clean
rm -rf build source/auto-kappa
sphinx-apidoc -f -e -o source/auto-kappa ../auto_kappa

#exit
#make clean

make html

