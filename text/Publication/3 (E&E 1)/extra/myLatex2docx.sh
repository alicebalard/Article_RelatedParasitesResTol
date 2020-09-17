#!/bin/sh
#$1=tex
#$2=bib

pandoc -s $1 --bibliography=$2 -o generated.docx
