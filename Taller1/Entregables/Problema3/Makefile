a:
	g++ Problema3_puntoa.cpp -o Sistema_Sol_Jupiter.out && ./Sistema_Sol_Jupiter.out
	cd plt && gnuplot SJ.p && cd ../
b:
	g++ Problema3_puntob.cpp -o rotado.out && ./rotado.out
	cd plt && gnuplot SJ_R.p && cd ../
c:
	g++ Problema3_puntoc.cpp -o troyano.out && ./troyano.out
	cd plt && gnuplot SJT.p && cd ../

d:
	g++ Problema3_puntod.cpp -o perturbado.out && ./perturbado.out
	cd plt && gnuplot SJT_p.p && cd ../

e:
	g++ Problema3_puntoe.cpp -o perturbado_periodos.out && ./perturbado_periodos.out
	cd plt && gnuplot SJT_p_periodos.p && cd ../
clean:
	rm -f *.out *~ *# *.dat
	clear
