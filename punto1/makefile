x-P.png : graficas.py UpwindGodunov_step_5.dat
	python graficas.py

UpwindGodunov_step_5.dat : shocktube1.x
	./shocktube1.x
	
shocktube1.x : shocktube1.c
	gcc shocktube1.c -o shocktube1.x -lm

# cuando se corre make la primera vez sale -- make: *** [UpwindGodunov_step_5.dat] Error 136 -- hacer make una segunda vez genera el resto de archivos
