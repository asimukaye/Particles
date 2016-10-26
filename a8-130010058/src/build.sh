mkdir -p ../output
python a8-130010058.py

jupyter nbconvert --execute --to html a8-130010058.ipynb

mv  *.png ../output

mv a8-130010058.html ../output/
mv a8-130010058.pdf ../output/
