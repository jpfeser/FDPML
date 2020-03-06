MODULE init

	USE kinds
	CONTAINS
	
!	===========================================================================================
	
	SUBROUTINE set_defaults(mass_input, crystal_coordinates, asr, mp, qpoint, &
						nk1, nk2, nk3, file_input, qlist_file, slist_file, Llist_file, &
						Lpoint, spoint, tol, maxit, restart, iterpause, calc_TC, &
						scattered_energy, plot_K, plot_sig, plot_uinc, plot_uscat)

		USE kinds
		IMPLICIT NONE
	
		CHARACTER(len = 256) 			:: 	mass_file, domain_file, simulation_type, &
											qlist_file, slist_file, Llist_file, asr
		LOGICAL							::	plot_K, plot_sig, plot_uinc, plot_uscat, calc_TC, &
											scattered_energy, restart, mass_input, crystal_coordinates, &
											mp, file_input
		INTEGER							::	iterpause
		REAL (KIND = RP) 				:: 	tol
		INTEGER(KIND = 4) 				:: 	maxit, flag
		INTEGER							::	spoint, Lpoint
		INTEGER							::	qpoint, nk1, nk2, nk3
		
		
	
	
	!	----------------------------			 
	!	Default values
	!	
		mass_input = .false.
		crystal_coordinates = .false.
		asr = 'simple'
		mp = .false.
		qpoint = 1
		nk1 = 4
		nk2 = 4
		nk3 = 4
		file_input = .false.
		qlist_file = 'NA'
		slist_file = 'NA'
		Llist_file = 'NA'
		Lpoint = 1
		spoint = 1
		tol = 1e-6
		maxit = 10000
		restart = .false.
		iterpause = 1000
		calc_TC = .false.
		scattered_energy = .false.
		plot_K = .false.
		plot_sig = .false.
		plot_uinc = .false.
		plot_uscat = .false.
	!
	!	--------------------------
	
	END SUBROUTINE
	
!	=============================================================================================
END MODULE init
