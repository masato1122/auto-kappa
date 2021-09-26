
PYTHONPATH=$PYTHONPATH:/Users/masato/work/aiida_alamode
export PYTHONPATH

rm *.in
python test_alm.py
python test_anphon.py

