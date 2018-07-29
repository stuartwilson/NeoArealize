#!/bin/bash

xelatex methods_notes.tex
bibtex methods_notes.aux
xelatex methods_notes.tex
xelatex methods_notes.tex
