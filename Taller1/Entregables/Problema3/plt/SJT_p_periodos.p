set terminal png size 500,500
set grid
set output '../img/e.png'
set title 'Sol Jupiter Troyano-perturbado Periodos'
plot "../S-J-T_perturbated_periodos.dat" w l title 'Troyano'
