figura16.png : graficas.py salidaLF.dat salida1LF.dat salida2LF.dat salida3LF.dat salidaRK4.dat salida1RK4.dat salida2RK4.dat salida3RK4.dat
	python graficas.py

salidaLF.dat : Integrador1.x 
	./Integrador1.x

salida1LF.dat : Integrador2.x 
	./Integrador2.x 

salida2LF.dat : Integrador3.x 
	./Integrador3.x 

salida3LF.dat : Integrador4.x 
	./Integrador4.x 

salidaRK4.dat : Integrador1.x 
	./Integrador1.x

salida1RK4.dat : Integrador2.x 
	./Integrador2.x 
	
salida2RK4.dat : Integrador3.x 
	./Integrador3.x
	
salida3RK4.dat : Integrador4.x
	./Integrador4.x

Integrador1.x : Integrador1.c
	gcc Integrador1.c -o Integrador1.x -lm

Integrador2.x : Integrador2.c
	gcc Integrador2.c -o Integrador2.x -lm
	
Integrador3.x : Integrador3.c
	gcc Integrador3.c -o Integrador3.x -lm
	
Integrador4.x : Integrador4.c
	gcc Integrador4.c -o Integrador4.x -lm
