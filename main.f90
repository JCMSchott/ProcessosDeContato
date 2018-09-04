!##################################################
!	Programa principal
!##################################################

program main
	use mod_rndgen
	use geraRede, only: grafoRRN
	use dynamics !, only: contactProcess
	
	implicit none
	
	!##########################################
	!	Parametros internos ao probrama
	!##########################################	
	integer, parameter :: dp = kind(0.0d0)

	!##########################################
	!	Medidor de tempo computacional
	!##########################################	

	real(dp) :: start, finish
	

	!##########################################
	!	Substrato ou rede e seus
	!	parametros
	!##########################################	


	type(grafoRRN) :: rede
	
	
	integer :: ki, N

	!##########################################
	!	Semente para o gerador de numeros
	!	pseudo aleatorios.
	!##########################################	

	integer :: seed
	
	!##########################################
	!	Inicializacao dos parametros
	!	da rede e do gerador de numeros
	!	pseudo aleatorios.
	!##########################################	

	N = 10**4
	ki = 4
	seed = 999999999
	
	
	call cpu_time(start)
	
		call rede%inicia(N, ki)

		write(*,*) " "
		write(*,*) "Rede inicializada"
		write(*,*) " "
		write(*,*) "Numero de nos= ", rede%nodes, "grau por no= ", rede%ki
		write(*,*) " "
		
		call rede%liga(seed)

		write(*,*) " "
		write(*,*) "Rede gerada com sucesso"
		write(*,*) " "
		
		
		lambda = 3.0;	nInf0 = 1d0
		
		call condicaoInicial(rede, nInf0, seed)
		call contactProcess(rede, tMax)

	call cpu_time(finish)
	write(*,*) "Tempo de execucao= ", finish - start

end program
