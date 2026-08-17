./sphere 2>&1 | tee out.log

cd ./data

python3 plot_p.py
