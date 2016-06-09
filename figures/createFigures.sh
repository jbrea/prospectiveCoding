#!/bin/bash
math -script src/mathematica/createFigures.m
cd figures
for i in Fig*[!0-9].tex; do latex $i; dvips ${i/tex/dvi}; ps2pdf ${i/tex/ps}; done
for i in Fig*[0-9].tex; do latex $i; dvips ${i/tex/dvi}; ps2pdf ${i/tex/ps}; done
latexmk -c
cd ..
