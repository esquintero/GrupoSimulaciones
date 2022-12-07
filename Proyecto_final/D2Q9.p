set title 'Plot de Adv_Dif'

set pm3d map interpolate 3,3

set size ratio -1

set palette rgbformulae 7,5,15

set cbrange [1:3]

splot "D2Q9.dat"

pause -1
