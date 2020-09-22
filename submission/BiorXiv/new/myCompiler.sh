#!/bin/sh
filename=$1

#echo $1
#echo "${filename%.*}"

xelatex $1
biber "${filename%.*}"
xelatex $1
xelatex $1
google-chrome "${filename%.*}".pdf
