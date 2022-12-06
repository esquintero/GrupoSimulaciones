f(x)=m*x+b

fit f(x) "varianza.dat" via m, b

plot "varianza.dat" with points ,f(x) title "Line Fit"

pause -1
 
