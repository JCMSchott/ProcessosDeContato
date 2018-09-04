!include 'mod_rndgen_multiple.f90'

!##################################################
!	Modulo que gera a rede
!##################################################

module	geraRede
	use mod_rndgen
	implicit none
	
	type grafo
		integer :: nodes, edge
		integer, allocatable:: deg(:), listAdj(:), aux(:), matriz(:,:)
		
		contains
			procedure :: iniciaGrafo
	end type
	
	type,extends(grafo) :: grafoRRN
		integer :: ki
	
		contains
			procedure ::	inicia => iniciaGrafoRRN
			procedure ::	liga => ligaRRN
	end type
	
	!###############################################################
	!		Subrotina que inicializa o grafo
	!###############################################################
	
	contains
	
		subroutine iniciaGrafo(this, N)
			integer, intent(in) :: N
			class(grafo) :: this
		
			this%nodes = N;	this%edge = 0
		
			if(allocated(this%deg)) deallocate(this%deg)
				allocate(this%deg(N)); this%deg = 0
			
			if(allocated(this%aux)) deallocate(this%aux)
				allocate(this%aux(N))
		end subroutine
		
		

		!###############################################################
		!		Setter pra inicializar junto ki.
		!		No inicializador padrao nao eh possivel.
		!###############################################################

		subroutine iniciaGrafoRRN(this, N, ki)
			class(grafoRRN) :: this
			integer, intent(in) :: ki, N
			integer :: i
			
			call iniciaGrafo(this,N)
			
			this%ki = ki
			
			this%aux(1) = 1
			
			do i = 2, this%nodes
				this%aux(i) = this%aux(i-1) + this%ki
			enddo
			
			if(allocated(this%listAdj)) deallocate(this%listAdj)
				allocate(this%listAdj(this%ki * this%nodes))
			
			if(allocated(this%matriz)) deallocate(this%matriz)
				allocate(this%matriz(this%nodes * this%ki/2, 2))
!				this%matriz = 0
		end subroutine



		!###############################################################
		!		Subrotina que conecta os nos do grafo RRN
		!###############################################################
	
	
		subroutine ligaRRN(this, seed)
			class(grafoRRN) :: this
			integer, allocatable :: lC(:)
			integer :: i, j, k
			integer :: seed
			integer :: cand1, cand2, lastC, nStubs
			logical :: floppou
			type(rndgen) :: gen
			
			open(10,file="conexoes.csv",status='unknown')
			write(10,*) "target",",", "source"
			call gen%init(seed)
			
			if(allocated(lC)) deallocate(lC)
				allocate(lC(this%ki * this%nodes))
								
			k = 1
			do i = 1, this%nodes
				do j = 1, this%ki
					lC(k) = i
					k = k + 1
				enddo
			enddo
			
			lastC = size(lC)
			nStubs = lastC/2
			write(*,*) "numero de stubs disponiveis", nStubs
			
do cand1 = 1, this%nodes
Principal:	do while(this%deg(cand1) < this%ki)
						cand2 = gen%int(1, lastC)
						if(this%deg(lC(cand2)) == this%ki)then
							cycle Principal
						endif
						if(lC(cand2) == cand1)then
							cycle Principal
						else
							floppou = .false.
							do k = this%aux(cand1), this%aux(cand1) + this%deg(cand1) - 1
								if(this%listAdj(k) == lC(cand2)) cycle Principal
							enddo
						endif
												
						this%listAdj(this%aux(cand1) + this%deg(cand1)) = lC(cand2)
						this%listAdj(this%aux(lC(cand2)) + this%deg(lC(cand2))) = cand1
						write(10, *) cand1, ",", lC(cand2)
						this%deg(cand1) = this%deg(cand1) + 1
						this%deg(lC(cand2)) = this%deg(lC(cand2)) + 1
						
						this%edge = this%edge + 1
						
						if(cand2 == lastC)then
							lastC = lastC - 1
						else
							lC(cand2) = lC(lastC); lastC = lastC - 1
						endif
			enddo Principal
enddo
			write(*,*) "Sobraram ", nStubs - this%edge," stubs"

			close(10)
		end subroutine
		
end module
