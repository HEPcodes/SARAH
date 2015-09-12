#!/bin/bash 
pdflatex Observable_PengV2L2u.tex 
cd Diagrams 
find . -name "*.mp" -exec mpost {} \; 
cd .. 
pdflatex  Observable_PengV2L2u.tex 
echo "" 
echo "PDF for Observable PengV2L2u finished" 
echo "Thanks for using SARAH" 
echo "" 
