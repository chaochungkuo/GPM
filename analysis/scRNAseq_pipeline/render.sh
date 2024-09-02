quarto render preprocessing/1_quality_control.ipynb --execute --to html && jupyter nbconvert --clear-output --inplace preprocessing/1_quality_control.ipynb
quarto render preprocessing/2_normalization.ipynb --execute --to html && jupyter nbconvert --clear-output --inplace preprocessing/2_normalization.ipynb
quarto render annotation/3_annotation.ipynb --execute --to html && jupyter nbconvert --clear-output --inplace annotation/3_annotation.ipynb

