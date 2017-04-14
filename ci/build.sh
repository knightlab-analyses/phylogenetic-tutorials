# build R files
input="r_requirements.txt"
while IFS= read -r var
do
  wget $var
  b=$(basename $var)
  R CMD INSTALL $b
done < "$input"
