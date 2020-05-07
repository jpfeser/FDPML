MODULE init

	USE kinds
	CONTAINS
	
!	===========================================================================================
	
	SUBROUTINE set_defaults(mass_input, crystal_coordinates, asr, qpoint, &
						file_input, qlist_file, slist_file, Llist_file, &
						Lpoint, spoint, tol, maxit, restart, iterpause, calc_TC, &
						scattered_energy, scattering_Xsec, plot_K, plot_sig, plot_uinc, plot_uscat, q_from_file, &
						q_file, expense_estimate)

		USE kinds
		IMPLICIT NONE
	
		CHARACTER(len = 256) 			:: 	mass_file, domain_file, simulation_type, &
											qlist_file, slist_file, Llist_file, asr
		LOGICAL							::	plot_K, plot_sig, plot_uinc, plot_uscat, calc_TC, &
											scattered_energy, restart, mass_input, crystal_coordinates, &
											mp, file_input, scattering_Xsec
		INTEGER							::	iterpause
		REAL (KIND = RP) 				:: 	tol
		INTEGER(KIND = 4) 				:: 	maxit, flag
		INTEGER							::	spoint, Lpoint
		INTEGER							::	qpoint, nk1, nk2, nk3
		LOGICAL							::	q_from_file, expense_estimate
		CHARACTER(len = 256)			::	q_file
		
		
	
	
	!	----------------------------			 
	!	Default values
	!	
		mass_input = .false.
		crystal_coordinates = .false.
		asr = 'simple'
		qpoint = 1
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
		q_from_file = .false.
		q_file = 'NA'
		scattering_Xsec = .false.
		expense_estimate = .false.
	!
	!	--------------------------
	
	END SUBROUTINE
	
!	=============================================================================================
END MODULE init
