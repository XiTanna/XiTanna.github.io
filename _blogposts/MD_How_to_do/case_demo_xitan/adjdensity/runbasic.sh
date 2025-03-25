

T=298
mdpname=grompp.mdp
groname=min.gro

export PATH=/share/home/xitan_n/soft/manycode/tuneDensity2020/changeMolNum:$PATH

## 2020
changeMolNM.sh   -n  4 -N 3552 -b 5000 -l 0.995 -u 1.005 -L 15.340 -U 17.340 -D   1048.58    -e 0.9 -f $mdpname -c $groname -p topol.top   



