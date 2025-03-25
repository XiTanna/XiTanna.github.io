lx=3.95360
ly=3.42240
lz=35.52
packmol<mixture.inp
gmx editconf -f mixture.pdb -o mixture.gro

sed -i "$ s/.*/$lx $ly $lz/" mixture.gro
rm -rf  ''#*.*.*#''