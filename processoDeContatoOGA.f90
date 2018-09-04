!########################################################################################################
!	Modulo que faz a dinamica
!########################################################################################################

module dynamics
	use geraRede
	use mod_rndgen
	implicit none
	
	!################################################################################################
	!	Parametros internos ao programa
	!################################################################################################

	integer, private, parameter :: dp = kind(0.0d0)
	
	!################################################################################################
	!	Listas dinamicas
	!################################################################################################
	
	integer, allocatable :: infList(:), susList(:), infSusList(:), sigma(:)


	!################################################################################################
	!	Variaveis dinamicas globais
	!################################################################################################

	
	real(dp), private, parameter :: mu = 1.0_dp
	
	real(dp) :: lambda
	

	!################################################################################################
	!	Para inicializar
	!################################################################################################
	
	real(dp) :: nInf0
	
	!################################################################################################
	!	Variaveis dinamicas
	!################################################################################################
	
	integer :: nInf, nSus, nInfSus

	
	!################################################################################################
	!	Taxas de eventos
	!################################################################################################

	real(dp) :: rateTotal, lambdaTotal, muTotal

	
	!################################################################################################
	!	Probabilidades associadas aos eventos independentes
	!################################################################################################
	
	real(dp) :: m, l
	

	!################################################################################################
	!	Variaveis temporais
	!################################################################################################

	real(dp) :: t, dt
	
	integer	:: t_pos, tMax


	!################################################################################################
	!	Variaveis epidemicas
	!################################################################################################

	
	real(dp), allocatable :: rho_medio(:), rho_medioQS(:)


	!################################################################################################
	!	gerador de numeros aleatorios
	!################################################################################################
	
	type(rndgen) :: gen

	
	!################################################################################################
	!	variavel aleatoria
	!################################################################################################
	
	real(dp) :: prob
	
	
	!################################################################################################
	!	Variaveis auxiliares
	!################################################################################################
	
	integer :: i, j, k




	contains
	
	
		!#######################################################################################
		!		Condicao inicial da dinamica
		!#######################################################################################
		
		
		subroutine condicaoInicial(this, nInf0, seed)
		
			!###############################################################################
			!	Variaveis auxiliares
			!###############################################################################
			
			integer :: lastCand
			
			integer, intent(in) :: seed
			
			!###############################################################################
			!	Variaveis inicializacao da dinamica
			!###############################################################################
			
			real(dp), intent(in) :: nInf0
			
			
			!###############################################################################
			!	A rede, ou substrato
			!###############################################################################

			class(grafo), intent(in) :: this

			!###############################################################################
			!	Lista aulixiar
			!###############################################################################
	
			integer, allocatable :: listAux(:)

			
		!#######################################################################################
		!			Situacao inicial da dinamica
		!#######################################################################################


			call gen%init(seed)
			
			nInf = int(nInf0 * this%nodes) 							! Quantidade de nos inicialmente infectados
			nSus = this%nodes - nInf							! Quantidade de nos inicialmente suscetiveis
													
			
			if(allocated(infList)) deallocate(infList)
				allocate(infList(this%nodes))
				

			if(allocated(sigma)) deallocate(sigma)
				allocate(sigma(this%nodes))
				sigma = 0								! A priori, todo mundo suscetivel
				
			if(nInf0 < 1)then								! Dava a porcentagem inicial de nos infectados										
				if(allocated(listAux)) deallocate(listAux)				! devemos inicializar a lista de nos inicialmente infectados
					allocate(listAux(this%nodes))

					do i = 1, this%nodes
						listAux(i) = i
					enddo
					
					k = 1;	lastCand = this%nodes
					do while(k <= nInf)

						i = gen%int(1, lastCand)				! lastCand pertence a geraRede
						infList(k) = listAux(i)
						sigma(listAux(i)) = 1

						listAux(i) = listAux(lastCand)
						lastCand = lastCand - 1

						k = k + 1
					enddo
					deallocate(listAux)
			else
				do i = 1, this%nodes
					infList(i) = i
					sigma(i) = 1
				enddo
			endif

!			lambdaTotal = 0.0_dp								! Aqui entram inclusive os processos
!			do i = 1, nInf									! fantasmas da epidemia, ou seja, 
			
!				lambdaTotal = nInf * lambda							
													! inclusive a tentativa de infectar
!			enddo										! quem ja esta infectado.
													
													! Se nao implementassemos o OGA, teriamos
													! que fazer:
													! lambdaTotal = lambdaTotal + &
													! 1.0_dp * kp(i) * lambda/this%deg(infList(i)),
													! onde kp seria uma lista contendo o numero
													! de vizinhos suscetiveis de infList.
													! Se infList(i) eh curado, ele eh retirado
													! de kp e de infList.
													! Mas no OGA isto fica mais simples.

!			muTotal = nInf * mu								! Taxa de cura inicial total			
!			rateTotal = muTotal + lambdaTotal						! Aqui calculamos inicialmente a taxa
													! total de eventos.
													! Ja atualizo essas taxas na proxima subrotina	
																																							
			!###############################################################################
			!	Finalizada aqui a situacao inicial, vamos partir para a dinamica
			!###############################################################################
			

			t = 0.0_dp; m = mu/(mu + lambda)
		
		end subroutine
		
		
		
		
		!#######################################################################################	
		!		Aqui segue a dinamica
		!#######################################################################################
		
		
		subroutine contactProcess(this, tMax)



			!###############################################################################
			!	A rede, ou substrato
			!###############################################################################

			class(grafo), intent(in) :: this
			
			real(dp), intent(in) :: tMax
			
			
loopDinamico:		do while(t <= tMax)

				muTotal = nInf * mu
				lambdaTotal = nInf * lambda
				rateTotal = muTotal + lambdaTotal				

				dt = -log((1.0_dp - gen%rnd())/rateTotal)/rateTotal
				t = t + dt
				
				prob = gen%rnd()
				if(prob <= m)then							! Se o processo da vez eh cura.
					i = gen%int(1, nInf)
					
					sigma(infList(i)) = 0
					infList(i) = infList(nInf)
					nInf = nInf - 1		
				else									! Se o processo da vez eh infeccao.
					i = gen%int(this%aux(i), this%aux(i) + this%deg(i) - 1)		! Posicao na lista de adjacencia correspondente aos vizinhos
													! de um sitio
					if(sigma(this%listAdj(i)) == 0)then								
						sigma(this%listAdj(i)) = 1				! Muda o estado do sitio.
						nInf = nInf + 1						! Numero de infectados aumenta.
						infList(nInf) = this%listAdj(i)				! Coloco esse novo infectado na primeira posicao disponivel
													! na lista de infectados.
					endif
				endif									


				
				if(nInf == 0) exit loopDinamico						! Porque o estado absorvente foi alcancado

			enddo loopDinamico			
						
			
			
			
		end subroutine
	
end module
