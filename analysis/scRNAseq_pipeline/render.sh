#!/bin/bash

EXECUTE_FLAG="--execute"
# EXECUTE_FLAG=""

# Preprocessing
cd preprocessing
pixi run quarto render 1_quality_control.ipynb $EXECUTE_FLAG --to html
pixi run jupyter nbconvert --clear-output --inplace 1_quality_control.ipynb
pixi run quarto render 2_normalization.ipynb $EXECUTE_FLAG --to html
pixi run jupyter nbconvert --clear-output --inplace 2_normalization.ipynb
cd ../

# Annotation
cd annotation
pixi run quarto render 3_annotation.ipynb $EXECUTE_FLAG --to html
pixi run jupyter nbconvert --clear-output --inplace 3_annotation.ipynb
cd ../