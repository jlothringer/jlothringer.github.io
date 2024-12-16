#!/bin/bash
source ~/.bash_profile
conda activate research
python jwst_target_page.py

cat Dashboard_Block025.html Dashboard_Block05.txt > tmp.txt
cat tmp.txt Dashboard_Block075.html > tmp2.txt
cat tmp2.txt Dashboard_Block2.html > tmp3.txt
cat tmp3.txt Dashboard_Block3.html > Dashboard.html

cp Dashboard.html ../
cp *png ../
cp *gif ../
cp JWST_ExoDashboard_Table*csv ../

rm tmp.txt
rm tmp2.txt
rm tmp3.txt

#git commit -a -m "update"
#git push origin main
