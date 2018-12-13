!  Tripe_decomposition_analysis.f90 

!****************************************************************************
!
!  PROGRAM: Tripe_decomposition_analysis
!
!  Version 2
!  by Stefan Felder - August 2009 - 2018
!
! in this version, the triple decomposition analyses are conducted with a factor! please note, that the mean calculations (i.e. low) are 
! standing for the calculations of correlations of slow+fast velocity fluctuations
!
!****************************************************************************

program Tripe_decomposition_analysis

implicit none

	!Parameters

	INTEGER, PARAMETER :: DBL = SELECTED_REAL_KIND(p=15)	!Precision of the Real(kind=DBL) 
															!variables
	!General variables
	!parameters to be loaded from file
	REAL :: delta_x
	INTEGER :: frequency
	INTEGER :: Correlation_steps
	INTEGER :: sample_duration
	INTEGER :: Number_devices
	INTEGER :: segment_numbers
	REAL :: PDF_V_segment_size
	INTEGER :: max_no_bubbles
	REAL :: threshold_C
	REAL :: PDF_Court_length_segment_size
	REAL :: PDF_Court_time_segment_size
	REAL :: PDF_Court_min_length
	REAL :: PDF_Court_max_length
	REAL :: PDF_Court_min_time
	REAL :: PDF_Court_max_time
	INTEGER :: cluster_or_not
	REAL :: cluster_factor
	REAL :: cluster_factor_1
	REAL :: lower_frequency
	REAL :: upper_frequency ! threshold in frequency)
	INTEGER :: FFT_number  ! integer for FFT analysis (2^n)
	INTEGER :: Correlation_steps_mid
	INTEGER :: Correlation_steps_low

	REAL(kind=DBL), ALLOCATABLE, DIMENSION(:,:) :: Data_array
	CHARACTER(len=9) :: name_parameters			!name of file with parameters to be opened
	CHARACTER(len=9) :: name					!name of file with positions to be opened
	CHARACTER(len=6) :: name_binary
	CHARACTER(len=8) :: location_binary
	CHARACTER(len=4) :: position
	CHARACTER(len=26) :: position_binary				
	CHARACTER(len=4) :: ending = '.txt'
	CHARACTER(len=1) :: sample = '_'
	INTEGER :: number_rows = 0					!number of rows to read
	INTEGER :: status = 0							!I/O status	
	REAL(kind=DBL) :: value	= 0
	INTEGER :: number_positions	= 0
	CHARACTER(len=26), ALLOCATABLE, DIMENSION(:) :: Position_array_binary
	CHARACTER(len=4), ALLOCATABLE, DIMENSION(:) :: Position_array
	INTEGER :: number_positions_binary  =0 
	CHARACTER(len=4) :: location					
	INTEGER :: number_columns = 0					
	INTEGER :: i,j,k

	!data arrays for annalyses

	! Variables for correlation-analyses
	REAL(kind=DBL), ALLOCATABLE, DIMENSION(:,:) :: Autocorrelation_raw			!results_file
	REAL(kind=DBL), ALLOCATABLE, DIMENSION(:,:) :: Crosscorrelation_raw		!results_file
	REAL(kind=DBL), ALLOCATABLE, DIMENSION(:,:) :: Autocorrelation_high			!results_file
	REAL(kind=DBL), ALLOCATABLE, DIMENSION(:,:) :: Crosscorrelation_high		!results_file
	REAL(kind=DBL), ALLOCATABLE, DIMENSION(:,:) :: Autocorrelation_mid			!results_file
	REAL(kind=DBL), ALLOCATABLE, DIMENSION(:,:) :: Crosscorrelation_mid		!results_file
	REAL(kind=DBL), ALLOCATABLE, DIMENSION(:,:) :: Autocorrelation_low			!results_file
	REAL(kind=DBL), ALLOCATABLE, DIMENSION(:,:) :: Crosscorrelation_low		!results_file

! for cross-products
	REAL(kind=DBL), ALLOCATABLE, DIMENSION(:,:) :: Autocorrelation_fast_low
	REAL(kind=DBL), ALLOCATABLE, DIMENSION(:,:) :: Crosscorrelation_fast_low
	REAL(kind=DBL), ALLOCATABLE, DIMENSION(:,:) :: Autocorrelation_low_fast
	REAL(kind=DBL), ALLOCATABLE, DIMENSION(:,:) :: Crosscorrelation_low_fast


	!Variables for Propability density distribution of Voltage signal for both probe tips
	REAL(kind=DBL) :: PDF_segments	= 0							!Number of segments for PDF
	REAL(kind=DBL), ALLOCATABLE, DIMENSION (: , :) :: PDF_V
	REAL(kind=DBL), ALLOCATABLE, DIMENSION(:) :: threshold			!results_file


	!Variables for Correlation-scale analyses
	REAL(kind=DBL), DIMENSION(15) :: CorScale_results_raw				!summary file with Correlation properties (raw)
	REAL(kind=DBL), DIMENSION(18) :: CorScale_results_high				!summary file with Correlation properties (high)
	REAL(kind=DBL), DIMENSION(18) :: CorScale_results_mid				!summary file with Correlation properties (mid)
	REAL(kind=DBL), DIMENSION(16) :: CorScale_results_low				!summary file with Correlation properties (mid)
	REAL(kind=DBL) :: velocity = 3									!default value

	!Variables for basic property analyses (air-water-interfaces, void fraction, frequency)
	INTEGER, ALLOCATABLE, DIMENSION(:,:) :: interfaces !interfaces
	REAL(kind=DBL), ALLOCATABLE, DIMENSION(:,:) :: basic_results	!C and F

	!Variables for chord time and size distributions
	REAL(kind=DBL) :: PDF_number_chord = 0
	REAL(kind=DBL) :: PDF_number_time = 0

	!variables for low and high pass filtering
	INTEGER :: Int_spectral_analys		! Integer for spectral analyses
	REAL(kind=DBL) :: sample_step		! time step between data points
	REAL(kind=DBL) :: total_FFT_duration
	REAL(kind=DBL) :: frequency_FFT
	INTEGER :: Cutoff_time_low_lower
	INTEGER :: Cutoff_time_low_upper
	INTEGER :: Cutoff_time_mid_lower
	INTEGER :: Cutoff_time_mid_upper
	INTEGER :: Cutoff_time_high_upper 
	INTEGER :: Cutoff_time_high_lower
	REAL(kind=DBL) :: helper_time
	REAL(kind=DBL) :: last_point_2
	REAL(kind=DBL) :: last_point_3
	REAL(kind=DBL) :: Total
	REAL(kind=DBL) :: Mean(3)
	REAL(kind=DBL), ALLOCATABLE, DIMENSION(:,:) :: FFT_array
		
	REAL(kind=DBL), ALLOCATABLE, DIMENSION(:,:) :: Filtered_array_low !(0-0.33Hz)
	REAL(kind=DBL), ALLOCATABLE, DIMENSION(:,:) :: Filtered_array_low_out !(0-0.33Hz)
	REAL(kind=DBL), ALLOCATABLE, DIMENSION(:,:) :: Filtered_array_mid !(0.33-10Hz)
	REAL(kind=DBL), ALLOCATABLE, DIMENSION(:,:) :: Filtered_array_high !(10-10000Hz)


	! data for correlation analyses of cut off data (3 seconds at each end cut off, e.g. 900000points- 2*60000 = 780000)

	INTEGER :: number_rows_cut = 0	
	REAL(kind=DBL) :: Auto_fac
	REAL(kind=DBL) :: Cross_fac
	REAL(kind=DBL) :: Factors_crossproduct(3)

	REAL(kind=DBL), ALLOCATABLE, DIMENSION(:,:) :: Data_array_cut  !(cut_array_for_correlation_analyses (3 sec at both ends cut off)
	REAL(kind=DBL), ALLOCATABLE, DIMENSION(:,:) :: Filtered_array_low_out_cut !(cut_array_for_correlation_analyses (3 sec at both ends cut off)
	REAL(kind=DBL), ALLOCATABLE, DIMENSION(:,:) :: Filtered_array_mid_cut    !(cut_array_for_correlation_analyses (3 sec at both ends cut off)
	REAL(kind=DBL), ALLOCATABLE, DIMENSION(:,:) :: Filtered_array_high_cut	!(cut_array_for_correlation_analyses (3 sec at both ends cut off)


	!Open file to get parameters from parameter file

	!Get the filename of the list of parameters for the analyses
	WRITE (*,*) 'Please input file name with parameters:'
	READ (*,*) name_parameters
	WRITE (*,501) name_parameters
	501 format (' ', 'The filename is: ', A '.txt')

	OPEN (UNIT = 21, FILE = name_parameters // ending, STATUS = 'OLD', ACTION ='READ', IOSTAT=status)
IF (status == 0) THEN

	READ (21, 100) delta_x
	100 FORMAT (//, 1F5.2)

	READ(21, 101) frequency
	101 FORMAT (/, 1I7)

	READ(21, 102) Correlation_steps
	102 FORMAT (/, 1I5)

	READ(21, 103) sample_duration
	103 FORMAT (//, 1I5)

	READ(21, 104) Number_devices
	104 FORMAT (//, 1I1)

	READ(21, 105) segment_numbers
	105 FORMAT (/, 1I4)

	READ(21, 106) PDF_V_segment_size
	106 FORMAT (//, 1F5.2)
			
	READ(21, 107) max_no_bubbles
	107 FORMAT (/, 1I6)

	READ(21, 108) threshold_C
	108 FORMAT (//, 1F5.2)

	READ(21, 109) PDF_Court_length_segment_size
	109 FORMAT (/, 1F5.2)

	READ(21, 110) PDF_Court_time_segment_size
	110 FORMAT (/, 1F5.2)

	READ(21, 111) PDF_Court_min_length
	111 FORMAT (/, 1F5.2)

	READ(21, 112) PDF_Court_max_length
	112 FORMAT (/, 1F5.2)

	READ(21, 113) PDF_Court_min_time
	113 FORMAT (/, 1F5.2)

	READ(21, 114) PDF_Court_max_time
	114 FORMAT (/, 1F5.2)

	READ(21, 115) 	cluster_or_not
	115 FORMAT (/, 1I1)

	READ(21, 116) cluster_factor
	116 FORMAT (//, 1F5.2)

	READ(21, 117) cluster_factor_1
	117 FORMAT (/, 1F5.2)

	READ(21, 118) lower_frequency
	118 FORMAT (/, 1F5.2)

	READ (21,119) upper_frequency
	119 FORMAT (/, 1F5.2)

	READ (21,120) FFT_number
	120 FORMAT (/, 1I8)

	READ (21,121) Correlation_steps_mid
	121 FORMAT (/, 1I5)

	READ (21,122) Correlation_steps_low
	122 FORMAT (/, 1I5)


END IF
    CLOSE ( UNIT=21)


	!Allocate arrays for results with parameters
	ALLOCATE (basic_results(Number_devices*2, 1), STAT =status)


	!Write Headings to result-files. Note: This needs to be done before the main loops-starts because the 
	!heading should just appear once in the sumary result files (Tu_V, basics and chordlength/-times
	OPEN (UNIT = 50, FILE =  'Tu_V_summary_raw' // ending, STATUS = 'REPLACE', &
		 ACTION ='WRITE', IOSTAT=status)
		WRITE (50, 150) 'Position', ' T0.5 ', ' T(Rxx=0) ', ' T(Min(Rxx)) ', ' Min(Rxx) ', ' Txx ', ' (Rxz)max ' , ' T(Rxz)max ', ' T(Rxz=0)&
			' , ' T(Min(Rxz)) ', ' Min(Rxz) ', ' T(0.5(Rxz)max) ' , ' Txz ' , ' Tu ' , ' V ', '(Tu*V)'
		150 FORMAT (A9, T11, A5, T21, A9, T31, A12, T43, A9, T52, A4, T62, A9, T72, A10, T83, A9, T94,&
			   A12, T106, A9, T115, A15, T131, A4, T141, A3, T151, A2, T161, A6)

	OPEN (UNIT = 94, FILE =  'Tu_V_summary_high' // ending, STATUS = 'REPLACE', &
		 ACTION ='WRITE', IOSTAT=status)
		WRITE (94, 1050) 'Position', ' T0.5 ', ' T(Rxx=0) ', ' T(Min(Rxx)) ', ' Min(Rxx) ', ' Txx ', ' (Rxz)max ' , ' T(Rxz)max ', ' T(Rxz=0)&
			' , ' T(Min(Rxz)) ', ' Min(Rxz) ', ' T(0.5(Rxz)max) ' , ' Txz ' , ' Tu ' , ' V ', '(Tu*V)', ' Rxx_max ', 'Rxx_fact', 'Rxz_fact'
		1050 FORMAT (A9, T11, A5, T21, A9, T31, A12, T43, A9, T52, A4, T62, A9, T72, A10, T83, A9, T94,&
			   A12, T106, A9, T115, A15, T131, A4, T141, A3, T151, A2, T161, A6, T172, A8, T182, A8, T192, A8)

	OPEN (UNIT = 95, FILE =  'Tu_V_summary_mid' // ending, STATUS = 'REPLACE', &
		 ACTION ='WRITE', IOSTAT=status)
		WRITE (95, 1051) 'Position', ' T0.5 ', ' T(Rxx=0) ', ' T(Min(Rxx)) ', ' Min(Rxx) ', ' Txx ', ' (Rxz)max ' , ' T(Rxz)max ', ' T(Rxz=0)&
			' , ' T(Min(Rxz)) ', ' Min(Rxz) ', ' T(0.5(Rxz)max) ' , ' Txz ' , ' Tu ' , ' V ', '(Tu*V)', ' Rxx_max ', 'Rxx_fact', 'Rxz_fact'
		1051 FORMAT (A9, T11, A5, T21, A9, T31, A12, T43, A9, T52, A4, T62, A9, T72, A10, T83, A9, T94,&
			   A12, T106, A9, T115, A15, T131, A4, T141, A3, T151, A2, T161, A6, T172, A8, T182, A8, T192, A8)

	OPEN (UNIT = 96, FILE =  'Tu_V_summary_sum_fast_slow_fluctuations' // ending, STATUS = 'REPLACE', &
		 ACTION ='WRITE', IOSTAT=status)
		WRITE (96, 1052) 'Position', ' T0.5 ', ' T(Rxx=0) ', ' T(Min(Rxx)) ', ' Min(Rxx) ', ' Txx ', ' (Rxz)max ' , ' T(Rxz)max ', ' T(Rxz=0)&
			' , ' T(Min(Rxz)) ', ' Min(Rxz) ', ' T(0.5(Rxz)max) ' , ' Txz ' , ' Tu ' , ' V ', '(Tu*V)', ' Rxx_max '
		1052 FORMAT (A9, T11, A5, T21, A9, T31, A12, T43, A9, T52, A4, T62, A9, T72, A10, T83, A9, T94,&
			   A12, T106, A9, T115, A15, T131, A4, T141, A3, T151, A2, T161, A6, T172, A8)

	OPEN (UNIT = 52, FILE =  'basics_summary'  // ending , STATUS = 'REPLACE', & 
		 ACTION ='WRITE', IOSTAT=status)
	WRITE (52, 152) 'Position', 'C-leading', 'F-leading', 'C-trailing', 'F-trailing'
		152 FORMAT (A9, T11, A10, T22, A10, T33, A11, T45, A11, T57)

! factor for crossproduct correlations for slow and fast fluctuations
	OPEN (UNIT = 520, FILE =  'factors_crossproduct_fast_slow'  // ending , STATUS = 'REPLACE', & 
		 ACTION ='WRITE', IOSTAT=status)
	WRITE (520, 252) 'Position', 'Auto_factor', 'Cross_slow_fast', 'Cross_fast_slow'
		252 FORMAT (A9, T11, A11, T23, A15, T39, A15, T54)

	number_columns = Number_devices + 1 			!Columns to be allocated in the arrays 
													!depending on the number of devices used
	number_rows = frequency * sample_duration		!Rows to be allocated in the data array

	PDF_segments = REAL(10.0)/PDF_V_segment_size		!Number of segments for PDF (segments between 0 and 10V)
	
	PDF_number_chord = 1+(PDF_Court_max_length-PDF_Court_min_length)/PDF_Court_length_segment_size

	PDF_number_time = 1+(PDF_Court_max_time-PDF_Court_min_time)/PDF_Court_time_segment_size


	!get the list with the binary file list
	WRITE (*,*) 'Please input file name with the binary file names:'
	READ (*,*) name_binary
	WRITE (*,510) name_binary
	510 format (' ', 'The name with the binary file names is  ', A '.txt')

	!Get the filename of the list of positions for the analyses
	WRITE (*,*) 'Please input file name with locations:'
	READ (*,*) name
	WRITE (*,500) name
	500 format (' ', 'The filename is: ', A '.txt')

	OPEN (UNIT = 20, FILE = name_binary // ending, STATUS = 'OLD', ACTION ='READ', IOSTAT=status)
	IF (status == 0) THEN

	!open was ok. Read values to find out how many lines are in the file.
		DO 
			READ (20, *, IOSTAT=status) location_binary	!Get next value
			IF (status /=0) EXIT						!EXIT if not valid
			number_positions_binary = number_positions_binary + 1		!Valid: increase count
		END DO 
	END IF

	OPEN (UNIT = 22, FILE = name // ending, STATUS = 'OLD', ACTION ='READ', IOSTAT=status)
	IF (status == 0) THEN

	!open was ok. Read values to find out how many lines are in the file.
		DO 
			READ (22, *, IOSTAT=status) location	!Get next value
			IF (status /=0) EXIT						!EXIT if not valid
			number_positions = number_positions + 1		!Valid: increase count	
		END DO 
	END IF

	IF (number_positions_binary == number_positions) THEN

	WRITE (*,520) number_positions
	520 format (' ', 'The file contains: ', I ' vertical locations in the cross section')

		ALLOCATE ( Position_array_binary(number_positions), STAT=status)	!allocate memory
		ALLOCATE ( Position_array(number_positions), STAT=status)

		! if the allocation was successful, rewind the file and read in the data!
		allocate_ok: IF ( status == 0) THEN
		REWIND ( UNIT=20)						!Rewind file	
		READ (20,*) (Position_array_binary(i), i=1, number_positions) 

		CLOSE ( UNIT=20)! Close File 
		END IF allocate_ok	

		allocate_ok2:IF ( status == 0) THEN
		REWIND ( UNIT=22)						!Rewind file	
		READ (22,*) (Position_array(i), i=1, number_positions)
		CLOSE ( UNIT=22)! Close File 
		END IF allocate_ok2	
	CLOSE ( UNIT=22)

	ELSE 
	WRITE (*,*) 'The number of positions in the two input files is different - please restart the program!'

	END IF	


DO i = 1, number_positions      !start of main loop in program!!!!
		position = Position_array(i)
		position_binary = Position_array_binary(i)			!Get the filename and echo it back to the user
		WRITE (*,10000) i, position
		10000 format (' ', 'The position ' I ' is: ' , A)

	! Open the position file, and check for errors on open.
	OPEN (UNIT = 30, FILE = position_binary // '.dat', STATUS = 'OLD', ACTION ='READ', IOSTAT=status)
	IF (status == 0) THEN

		! ALLOCATE MEMORY		
		ALLOCATE ( Data_array(number_columns, number_rows), STAT=status)	!allocate memory
		Data_array = 0

		! if the allocation was successful, rewind the file and read in the data!
		allocate_ok3: IF ( status == 0) THEN

		READ (30,*) ((Data_array (j,k), j=2,number_columns), k=1,(number_rows))

		END IF allocate_ok3				
		END IF	 
		
	! Close File 
	CLOSE ( UNIT=30)

	!Calculate basic data using the full raw signal
	! allocata arrays for analyses
	ALLOCATE (threshold(Number_devices),STAT=status)
	ALLOCATE (interfaces(Number_devices*2, max_no_bubbles), STAT =status)


	!Calculate the Probability distribution function of the voltage signals and store to array PDF_V
	ALLOCATE (PDF_V(number_columns, Int(PDF_segments)), STAT=status)	!allocate memory for PDF_V
			PDF_V = 0
	CALL Probability_V (Data_array, number_columns, number_rows, PDF_V_segment_size, &
						PDF_segments, PDF_V, threshold, threshold_C)

	CALL Basic_properties (Data_array, number_columns, number_rows, threshold,&
						 Interfaces, basic_results, max_no_bubbles, frequency)

	!Write PDF_V distribution data to file: PDF_V-results
	OPEN (UNIT = 64, FILE = position // sample //'PDF_V'  // ending, STATUS = 'REPLACE', & 
		 ACTION ='WRITE', IOSTAT=status)
	WRITE (64,164) 'Bin-size', 'V-leading', 'V-trailing'
	WRITE (64, 165) ((PDF_V (j,k), j=1,3), k=1, INT(PDF_segments))
	164 FORMAT (A8, T13, A9, T25, A10)
	165 FORMAT (F10.5, T13, F10.7, T25, F10.7)

	DEALLOCATE (PDF_V, STAT = status)
				
	OPEN (UNIT = 66, FILE = position //sample // 'interfaces'  // ending , STATUS = 'REPLACE', & 
		 ACTION ='WRITE', IOSTAT=status)
	WRITE (66,166) 'W->A leading', 'A->W leading', 'W->A trailing', 'A->W trailing'
	WRITE (66, 167) ((interfaces (j,k), j=1,((number_columns-1)*2)), k=1, 10000)
	166 FORMAT (A12, T15, A12, T36, A13, T52, A13)
	167 FORMAT (I, T15, I, T36, I, T52, I)

	DEALLOCATE (interfaces, STAT =status)
	DEALLOCATE (threshold,STAT=status)

	!Write to file basic_summary
	WRITE (52, 153) position, ((basic_results (j,k), j=1,(number_columns-1)*2), k=1,1)
	153 FORMAT (A4, T9, F10.5, T20, F10.5, T31, F10.5, T44, F10.5, T55)

	
! Analyses of filtered signals for triple decomposition approach	
! Filtering of the raw data for tripledecomposition

	ALLOCATE (Filtered_array_low(number_columns, number_rows), STAT =status)
	ALLOCATE (Filtered_array_low_out(number_columns, number_rows), STAT =status)
	ALLOCATE (Filtered_array_mid(number_columns, number_rows), STAT =status)
	ALLOCATE (Filtered_array_high(number_columns, number_rows), STAT =status)
	Filtered_array_low = 0
	Filtered_array_low_out = 0
	Filtered_array_mid = 0
	Filtered_array_high = 0	


	helper_time = 0
Do j = 1, number_rows
		Filtered_array_low (1, j) = helper_time
		Filtered_array_low_out (1, j) = helper_time
		Filtered_array_mid (1, j) = helper_time
		Filtered_array_high (1, j) = helper_time
		helper_time = helper_time + 1.0/REAL(frequency) 
END DO
	
	Total = 0.0
	Mean = 0.0

	Do k =2, number_columns
		Do j = 1, number_rows
			Total = Total + Data_array (k,j)
		END DO
			Mean(k) = Total/REAL(number_rows)
			Total = 0.0
	END DO

	!	add zeroes so that data points=2^n before filterinG and subtract Mean values
	ALLOCATE (FFT_array	(number_columns, FFT_number), STAT=status)
	FFT_array = 0

	Do k = 2, number_columns
		Do j = 1, number_rows
			FFT_array (k,j) = Data_array (k,j) - Mean(k)
		END DO
	END DO

!	Do k = 2, number_columns
		Do j = number_rows + 1, FFT_number
			FFT_array (2,j) = 0
			FFT_array (3,j) = 0
		END DO
!	END DO

!    === Numerical filtering and Spectral analysis ===

	sample_step = 1/REAL(frequency)
	Int_spectral_analys = FFT_number/2
    total_FFT_duration = sample_step * REAL(FFT_number)
    frequency_FFT = 1./total_FFT_duration


	Cutoff_time_low_lower = INT((0.0 + frequency_FFT/2.)/frequency_FFT)+1
	Cutoff_time_low_upper = Int((lower_frequency - frequency_FFT/2.)/frequency_FFT)+1
	Cutoff_time_mid_lower = Int((lower_frequency + frequency_FFT/2.)/frequency_FFT)+1
	Cutoff_time_mid_upper = Int((upper_frequency - frequency_FFT/2.)/frequency_FFT)+1
	Cutoff_time_high_lower = Int((upper_frequency + frequency_FFT/2.)/frequency_FFT)+1
	Cutoff_time_high_upper = Int((frequency/2 - frequency_FFT/2.)/frequency_FFT)+1

	CALL Filter (number_rows, number_columns, FFT_number, Int_spectral_analys, FFT_array, Cutoff_time_low_lower, &
				Cutoff_time_low_upper, Cutoff_time_mid_lower, Cutoff_time_mid_upper, Cutoff_time_high_lower, &
				Cutoff_time_high_upper,	Filtered_array_low, Filtered_array_mid, Filtered_array_high)

	Do k = 2, number_columns
		Do j = 1, number_rows
			Filtered_array_low_out (k,j) = Filtered_array_low (k,j) + Mean(k)
		END DO
	END DO

	DEALLOCATE (Filtered_array_low, STAT =status)
	DEALLOCATE (FFT_array, STAT =status)
	

	!Write  low filtered data to file: filter-results
	OPEN (UNIT = 6000, FILE = position // sample // 'filtered_low' //  ending, STATUS = 'REPLACE', & 
		 ACTION ='WRITE', IOSTAT=status)
	WRITE (6000, 1600) 'time', 'lead', 'trai'
	WRITE (6000, 1610) ((Filtered_array_low_out (j,k), j=1,3), k=1, number_rows)
	1600 FORMAT (A4, T10, A4, T22, A4)
	1610 FORMAT (F10.6,T10, F10.5, T22, F10.7)	
	CLOSE (UNIT=6000)

	!Write  mid filtered data to file: filter-results
	OPEN (UNIT = 6001, FILE = position // sample // 'filtered_mid' //  ending, STATUS = 'REPLACE', & 
		 ACTION ='WRITE', IOSTAT=status)
	WRITE (6001, 1601) 'time', 'lead', 'trai'
	WRITE (6001, 1611) ((Filtered_array_mid (j,k), j=1,3), k=1, number_rows)
	1601 FORMAT (A4, T10, A4, T22, A4)
	1611 FORMAT (F10.6,T10, F10.5, T22, F10.7)	
	CLOSE (UNIT=6001)

	!Write  high filtered data to file: filter-results
	OPEN (UNIT = 6002, FILE = position // sample // 'filtered_high' //  ending, STATUS = 'REPLACE', & 
		 ACTION ='WRITE', IOSTAT=status)
	WRITE (6002, 1602) 'time', 'lead', 'trai'
	WRITE (6002, 1612) ((Filtered_array_high (j,k), j=1,3), k=1, number_rows)
	1602 FORMAT (A4, T10, A4, T22, A4)
	1612 FORMAT (F10.6,T10, F10.5, T22, F10.7)	
	CLOSE (UNIT=6002)

	number_rows_cut = number_rows - (6 * 20000)


	!Cut_off of data for correlation analyses	
	ALLOCATE (Data_array_cut(number_columns, number_rows_cut), STAT =status)
	ALLOCATE (Filtered_array_low_out_cut(number_columns, number_rows_cut), STAT =status)
	ALLOCATE (Filtered_array_mid_cut(number_columns, number_rows_cut), STAT =status)
	ALLOCATE (Filtered_array_high_cut(number_columns, number_rows_cut), STAT =status)
	Data_array_cut = 0
	Filtered_array_low_out_cut = 0
	Filtered_array_mid_cut = 0
	Filtered_array_high_cut = 0

	Do k = 1, number_columns
		Do j = 1, number_rows_cut 
			Data_array_cut (k, j)  = Data_array (k, j + (60000))
			Filtered_array_low_out_cut (k, j) = Filtered_array_low_out (k, j + (60000))
			Filtered_array_mid_cut (k, j) = Filtered_array_mid (k, j + (60000))
			Filtered_array_high_cut (k, j) = Filtered_array_high (k, j + (60000))
		End Do
	END DO

		DEALLOCATE (Data_array, STAT = status)
		DEALLOCATE (Filtered_array_low_out, STAT =status)
		DEALLOCATE (Filtered_array_mid, STAT =status)
		DEALLOCATE (Filtered_array_high, STAT =status)

	! perform correlation analyses for raw data part
	!Allocate arrays for correlation analyses
	ALLOCATE (Autocorrelation_raw(2, Correlation_steps), STAT =status)
	ALLOCATE (Crosscorrelation_raw(2, Correlation_steps), STAT =status)	

	!use the data arry and do auto- and cross -correlation analyses
	CALL Correlation_raw(Data_array_cut, number_columns, number_rows_cut, Correlation_steps, & 
			         Autocorrelation_raw, Crosscorrelation_raw, frequency, segment_numbers)

	!calculate the characteristic scales, Tu and Interfacial Velocity
	CALL Correlation_scales_raw (Autocorrelation_raw, Crosscorrelation_raw, Correlation_steps, & 
			 CorScale_results_raw, frequency, delta_x, velocity) 

	!Write autocorrelation data to file: Autocorrelation-results
	OPEN (UNIT = 60, FILE = position // sample // 'autocor_raw' //  ending, STATUS = 'REPLACE', & 
		 ACTION ='WRITE', IOSTAT=status)
	WRITE (60, 160) 'Time', 'Rxx'
	WRITE (60, 161) ((Autocorrelation_raw (j,k), j=1,2), k=1, Correlation_steps)
	160 FORMAT (A4, T15, A3)
	161 FORMAT (F10.5,T15, F10.7)

	!Write crosscorrelation data to file: Crosscorrelation-results
	OPEN (UNIT = 62, FILE = position //sample //'crosscor_raw'  // ending, STATUS = 'REPLACE', & 
		 ACTION ='WRITE', IOSTAT=status)
	WRITE (62, 162) 'Time', 'Rxy'
	WRITE (62, 163) ((Crosscorrelation_raw (j,k), j=1,2), k=1, Correlation_steps)
	162 FORMAT (A4, T18, A3)
	163 FORMAT (F14.9, T16, F14.9)

	CLOSE (UNIT=60)
	CLOSE (UNIT=62)

	DEALLOCATE (Autocorrelation_raw, STAT =status)
	DEALLOCATE (Crosscorrelation_raw, STAT =status)	

	!Write data to file: Turbulence and velocity results
	WRITE (50, 151) Position, (CorScale_results_raw )	
	151 FORMAT (A4, T11, F9.6, T21, F9.6, T31, F9.6, T43, F9.6, T52, F9.6, T62, F9.6, T72, &
			   F9.6, T83, F9.6, T94, F9.6, T106, F9.6, T115, F9.6, T131, F9.6, T141, F9.6, T151, F9.6, T161, F9.6)

	DEALLOCATE (Data_array_cut, STAT = status)

	!--------------------------------------------------------------
	! perform correlation analyses for high fluctuating data part
	!Allocate arrays for correlation analyses
	ALLOCATE (Autocorrelation_high(2, Correlation_steps), STAT =status)
	ALLOCATE (Crosscorrelation_high(2, Correlation_steps), STAT =status)	

	!use the data arry and do auto- and cross -correlation analyses
	CALL Correlation_high(Filtered_array_high_cut, Filtered_array_mid_cut, number_columns, number_rows_cut, Correlation_steps, & 
			         Autocorrelation_high, Crosscorrelation_high, frequency, segment_numbers, Auto_fac, Cross_fac)

	velocity = 3

	!calculate the characteristic scales, Tu and Interfacial Velocity
	CALL Correlation_scales_high (Autocorrelation_high, Crosscorrelation_high, Correlation_steps, & 
			 CorScale_results_high, frequency, delta_x, velocity, Auto_fac, Cross_fac) 

	!Write autocorrelation data to file: Autocorrelation-results
	OPEN (UNIT = 90, FILE = position // sample // 'autocor_high' //  ending, STATUS = 'REPLACE', & 
		 ACTION ='WRITE', IOSTAT=status)
	WRITE (90, 190) 'Time', 'Rxx'
	WRITE (90, 191) ((Autocorrelation_high (j,k), j=1,2), k=1, Correlation_steps)
	190 FORMAT (A4, T15, A3)
	191 FORMAT (F10.5,T15, F10.7)

	!Write crosscorrelation data to file: Crosscorrelation-results
	OPEN (UNIT = 92, FILE = position //sample //'crosscor_high'  // ending, STATUS = 'REPLACE', & 
		 ACTION ='WRITE', IOSTAT=status)
	WRITE (92, 192) 'Time', 'Rxy'
	WRITE (92, 193) ((Crosscorrelation_high (j,k), j=1,2), k=1, Correlation_steps)
	192 FORMAT (A4, T18, A3)
	193 FORMAT (F14.9, T16, F14.9)

	CLOSE (UNIT=90)
	CLOSE (UNIT=92)

	!Write data to file: Turbulence and velocity results
	WRITE (94, 194) Position, (CorScale_results_high )	
	194 FORMAT (A4, T11, F9.6, T21, F9.6, T31, F9.6, T43, F9.6, T52, F9.6, T62, F9.6, T72, &
			   F9.6, T83, F9.6, T94, F9.6, T106, F9.6, T115, F9.6, T131, F9.6, T141, F9.6, &
			   T151, F9.6, T161, F9.6, T172, F9.6, T181, F9.6, T191, F9.6)

	Auto_fac = 0
	Cross_fac = 0

	!--------------------------------------------------------------
	! perform correlation analyses for slow fluctuating data part
	!Allocate arrays for correlation analyses

	ALLOCATE (Autocorrelation_mid(2, Correlation_steps), STAT =status)
	ALLOCATE (Crosscorrelation_mid(2, Correlation_steps), STAT =status)	

	!use the data arry and do auto- and cross -correlation analyses
	CALL Correlation_mid(Filtered_array_mid_cut, Filtered_array_high_cut, number_columns, number_rows_cut, Correlation_steps, & 
			         Autocorrelation_mid, Crosscorrelation_mid, frequency, segment_numbers, Auto_fac, Cross_fac)

	velocity = 3

	!calculate the characteristic scales, Tu and Interfacial Velocity
	CALL Correlation_scales_mid (Autocorrelation_mid, Crosscorrelation_mid, Correlation_steps, & 
			 CorScale_results_mid, frequency, delta_x, velocity, Auto_fac, Cross_fac) 

	!Write autocorrelation data to file: Autocorrelation-results
	OPEN (UNIT = 100, FILE = position // sample // 'autocor_mid' //  ending, STATUS = 'REPLACE', & 
		 ACTION ='WRITE', IOSTAT=status)
	WRITE (100, 290) 'Time', 'Rxx'
	WRITE (100, 291) ((Autocorrelation_mid (j,k), j=1,2), k=1, Correlation_steps)
	290 FORMAT (A4, T15, A3)
	291 FORMAT (F10.5,T15, F10.7)

	!Write crosscorrelation data to file: Crosscorrelation-results
	OPEN (UNIT = 102, FILE = position //sample //'crosscor_mid'  // ending, STATUS = 'REPLACE', & 
		 ACTION ='WRITE', IOSTAT=status)
	WRITE (102, 292) 'Time', 'Rxy'
	WRITE (102, 293) ((Crosscorrelation_mid (j,k), j=1,2), k=1, Correlation_steps)
	292 FORMAT (A4, T18, A3)
	293 FORMAT (F14.9, T16, F14.9)

	CLOSE (UNIT=100)
	CLOSE (UNIT=102)

	!Write data to file: Turbulence and velocity results
	WRITE (95, 294) Position, (CorScale_results_mid )	
	294 FORMAT (A4, T11, F9.6, T21, F9.6, T31, F9.6, T43, F9.6, T52, F9.6, T62, F9.6, T72, &
			   F9.6, T83, F9.6, T94, F9.6, T106, F9.6, T115, F9.6, T131, F9.6, T141, F9.6, &
			   T151, F9.6, T161, F9.6, T172, F9.6, T181, F9.6, T191, F9.6)


	!--------------------------------------------------------------

	!--------------------------------------------------------------
	! add slow and fast fluctuating correlation function and perform analyses 
	!Allocate arrays

	ALLOCATE (Autocorrelation_low(2, Correlation_steps), STAT =status)
	ALLOCATE (Crosscorrelation_low(2, Correlation_steps), STAT =status)	

	DO j = 1, Correlation_steps
		Autocorrelation_low (1, j) = Autocorrelation_high (1, j)
		Autocorrelation_low (2, j) = Autocorrelation_high (2, j) + Autocorrelation_mid (2, j)
		Crosscorrelation_low (1, j) = Crosscorrelation_high (1, j)
		Crosscorrelation_low (2, j) = Crosscorrelation_high (2, j) + Crosscorrelation_mid (2, j)
	END DO

	velocity = 3

	!calculate the characteristic scales, Tu and Interfacial Velocity
	CALL Correlation_scales_low (Autocorrelation_low, Crosscorrelation_low, Correlation_steps, & 
			 CorScale_results_low, frequency, delta_x, velocity) 

	!Write autocorrelation data to file: Autocorrelation-results
	OPEN (UNIT = 103, FILE = position // sample // 'autocor_sum_slow_fast_fluctuations' //  ending, STATUS = 'REPLACE', & 
		 ACTION ='WRITE', IOSTAT=status)
	WRITE (103, 390) 'Time', 'Rxx'
	WRITE (103, 391) ((Autocorrelation_low (j,k), j=1,2), k=1, Correlation_steps)
	390 FORMAT (A4, T15, A3)
	391 FORMAT (F10.5,T15, F10.7)

	!Write crosscorrelation data to file: Crosscorrelation-results
	OPEN (UNIT = 104, FILE = position //sample //'crosscor_sum_slow_fast_fluctuations'  // ending, STATUS = 'REPLACE', & 
		 ACTION ='WRITE', IOSTAT=status)
	WRITE (104, 396) 'Time', 'Rxy'
	WRITE (104, 397) ((Crosscorrelation_low (j,k), j=1,2), k=1, Correlation_steps)
	396 FORMAT (A4, T18, A3)
	397 FORMAT (F14.9, T16, F14.9)

	CLOSE (UNIT=104)
	CLOSE (UNIT=103)
	

	!Write data to file: Turbulence and velocity results
	WRITE (96, 295) Position, (CorScale_results_low )	
	295 FORMAT (A4, T11, F9.6, T21, F9.6, T31, F9.6, T43, F9.6, T52, F9.6, T62, F9.6, T72, &
			   F9.6, T83, F9.6, T94, F9.6, T106, F9.6, T115, F9.6, T131, F9.6, T141, F9.6, T151, F9.6, T161, F9.6, T172, F9.6)

	DEALLOCATE (Autocorrelation_low, STAT =status)
	DEALLOCATE (Crosscorrelation_low, STAT =status)	
	DEALLOCATE (Autocorrelation_mid, STAT =status)
	DEALLOCATE (Crosscorrelation_mid, STAT =status)	
	DEALLOCATE (Autocorrelation_high, STAT =status)
	DEALLOCATE (Crosscorrelation_high, STAT =status)

!-----------------------------------------------------------------	
	! calculate auto-and cross correlation function for crossproduct of fast and slow fluctuating components Rx'y'' and Rx''y'

	ALLOCATE (Autocorrelation_fast_low(2, Correlation_steps_low), STAT =status)
	ALLOCATE (Crosscorrelation_fast_low(2, Correlation_steps_low), STAT =status)	
	ALLOCATE (Autocorrelation_low_fast(2, Correlation_steps_low), STAT =status)
	ALLOCATE (Crosscorrelation_low_fast(2, Correlation_steps_low), STAT =status)	

	!use the data arry and do auto- and cross -correlation analyses
	CALL Correlation_low(Filtered_array_mid_cut, Filtered_array_high_cut, number_columns, number_rows_cut, Correlation_steps_low, & 
			         Autocorrelation_fast_low, Crosscorrelation_fast_low, Autocorrelation_low_fast, Crosscorrelation_low_fast, &
					 frequency, segment_numbers, Factors_crossproduct)


	!Write autocorrelation data to file: Autocorrelation-results
	OPEN (UNIT = 113, FILE = position // sample // 'autocor_crossproduct_fast_slow_fluctuations' //  ending, STATUS = 'REPLACE', & 
		 ACTION ='WRITE', IOSTAT=status)
	WRITE (113, 310) 'Time', 'Rxx'
	WRITE (113, 311) ((Autocorrelation_fast_low (j,k), j=1,2), k=1, Correlation_steps_low)
	310 FORMAT (A4, T15, A3)
	311 FORMAT (F10.5,T15, F10.7)


	OPEN (UNIT = 123, FILE = position // sample // 'autocor_crossproduct_slow_fast_fluctuations' //  ending, STATUS = 'REPLACE', & 
		 ACTION ='WRITE', IOSTAT=status)
	WRITE (123, 320) 'Time', 'Rxx'
	WRITE (123, 321) ((Autocorrelation_low_fast (j,k), j=1,2), k=1, Correlation_steps_low)
	320 FORMAT (A4, T15, A3)
	321 FORMAT (F10.5,T15, F10.7)

	!Write crosscorrelation data to file: Crosscorrelation-results
	OPEN (UNIT = 114, FILE = position //sample //'crosscor_crossproduct_fast_slow_fluctuations'  // ending, STATUS = 'REPLACE', & 
		 ACTION ='WRITE', IOSTAT=status)
	WRITE (114, 316) 'Time', 'Rxy'
	WRITE (114, 317) ((Crosscorrelation_fast_low (j,k), j=1,2), k=1, Correlation_steps_low)
	316 FORMAT (A4, T18, A3)
	317 FORMAT (F14.9, T16, F14.9)

	OPEN (UNIT = 124, FILE = position //sample //'crosscor_sum_crossproduct_slow_fast_fluctuations'  // ending, STATUS = 'REPLACE', & 
		 ACTION ='WRITE', IOSTAT=status)
	WRITE (124, 326) 'Time', 'Rxy'
	WRITE (124, 327) ((Crosscorrelation_low_fast (j,k), j=1,2), k=1, Correlation_steps_low)
	326 FORMAT (A4, T18, A3)
	327 FORMAT (F14.9, T16, F14.9)


	WRITE (520, 253) Position, (Factors_crossproduct)
		253 FORMAT (A4, T10, F9.6, T22, F9.6, T38, F9.6, T53)


	CLOSE (UNIT=113)
	CLOSE (UNIT=114)
	CLOSE (UNIT=123)
	CLOSE (UNIT=124)

	DEALLOCATE (Autocorrelation_fast_low, STAT =status)
	DEALLOCATE (Crosscorrelation_fast_low, STAT =status)	
	DEALLOCATE (Autocorrelation_low_fast, STAT =status)
	DEALLOCATE (Crosscorrelation_low_fast, STAT =status)	


	
	DEALLOCATE (Filtered_array_low_out_cut, STAT =status)
	DEALLOCATE (Filtered_array_mid_cut, STAT =status)
	DEALLOCATE (Filtered_array_high_cut, STAT =status)

!--------------------------------------------------------------

		WRITE (*,6000) position
		6000 format (' ', 'Finished calculation of position: ' A)

END DO				! END of main loop in program

	
	DEALLOCATE (basic_results, STAT =status)

	CLOSE (UNIT=50)
	CLOSE (UNIT=52)
	CLOSE (UNIT=54)
	CLOSE (UNIT=56)
	CLOSE (UNIT=64)
	CLOSE (UNIT=66)
	CLOSE (UNIT=94)
	CLOSE (UNIT=95)
	CLOSE (UNIT=96)
	CLOSE (UNIT=6000)
	CLOSE (UNIT=6001)
	CLOSE (UNIT=6002)
	CLOSE (UNIT=520)

	pause

end program dataanalysis

!______________________________________________________________________________________________________

SUBROUTINE Correlation_raw (Data_array_cut, number_columns, number_rows_cut, Correlation_steps, &  
			                Autocorrelation_raw, Crosscorrelation_raw, frequency, segment_numbers)
IMPLICIT NONE

	INTEGER, PARAMETER :: DBL = SELECTED_REAL_KIND(p=15)	!Precision of the Real variables
	INTEGER, INTENT (IN) :: number_columns
	INTEGER, INTENT (IN) :: number_rows_cut
	INTEGER, INTENT (IN) :: Correlation_steps
	INTEGER, INTENT (IN) :: frequency
	INTEGER, INTENT (IN) :: segment_numbers
	REAL(kind=DBL), INTENT (IN) :: Data_array_cut (number_columns, number_rows_cut)
	REAL(kind=DBL), INTENT (OUT) :: Autocorrelation_raw (2, Correlation_steps)
	REAL(kind=DBL), INTENT (OUT) :: Crosscorrelation_raw(2, Correlation_steps)
	INTEGER :: segment_rows							!number of rows in each segment
	INTEGER :: helper								!helper for faster correlations
	INTEGER :: boundary								!helper for faster correlation
	INTEGER :: i =0, j = 0, k=0							!integer for loop
	REAL(kind=DBL), ALLOCATABLE, DIMENSION  (:,:) :: Autocor_helper	
	REAL(kind=DBL), ALLOCATABLE, DIMENSION  (:,:) :: Crosscor_helper	
	REAL(kind=DBL) :: sumX, sumY, sumX2, sumY2, sumXY !Variables for correl
	INTEGER :: offset												 !integer for offsetting
	INTEGER :: status
	REAL(kind=DBL) :: average_value				!Helper to calculate the average values 												
	
	Autocorrelation_raw = 0
	Crosscorrelation_raw = 0
	segment_rows	=0
	helper = 0
	boundary = 0
	sumX = 0
	sumY = 0
	sumX2 = 0
	sumY2 = 0
	sumXY = 0
	offset = 0
	status = 0
	average_value =0


	segment_rows = number_rows_cut/segment_numbers
	boundary = segment_rows - Correlation_steps	


	!Performing Auto-correlation
	ALLOCATE (Autocor_helper(segment_numbers, Correlation_steps))
	Autocor_helper = 0			
	outer : DO i=1, segment_numbers
			helper = segment_rows * (i-1)
		middle : DO j=1, Correlation_steps
			inner : DO k = 1, boundary
					sumX = sumX + Data_array_cut(2 , (k + helper))
					sumX2 = sumX2 + Data_array_cut(2 , (k+ helper))* &
							Data_array_cut(2 , (k+ helper)) 
					sumY = sumY + Data_array_cut(2, (k+ helper+offset))
					sumY2 = sumY2 + Data_array_cut(2, (k+ helper+offset)) &
							*Data_array_cut(2, (k+ helper+offset))
					sumXY = sumXY + Data_array_cut(2,(k+ helper)) * &
							Data_array_cut(2,(k+ helper+offset))
			END DO inner
		
			Autocor_helper(i,j) = (sumXY - sumX * sumY/boundary) /(SQRT(sumX2-sumX*sumX/ &
								  boundary) * SQRT(sumY2-sumY*sumY/boundary))
			! Adjust variables
			offset=offset +1
			sumX = 0
			sumX2 = 0
			sumY = 0 
			sumY2 = 0 
			sumXY = 0 			

		END DO middle
		offset = 0		!zero offset helper
	END DO outer
	helper = 0

	!average the autocorrelation-segments and add time
	DO i = 1, Correlation_steps
			Autocorrelation_raw (1, i) = real((i-1))/frequency
			
			DO j = 1, segment_numbers
					average_value = average_value + Autocor_helper (j , i)
			END DO 
				
				Autocorrelation_raw (2, i) = average_value/segment_numbers
				average_value = 0
				
	END DO

	

	WRITE (*,*) 'Finished calculation of autocorrelation_raw! '

!		OPEN (UNIT = 2002, FILE = 'auto.txt', STATUS = 'REPLACE', ACTION ='WRITE', IOSTAT=status)
!		WRITE (2002, 2003) ((Autocor_helper (i,j), i=1,segment_numbers), j=1, Correlation_steps)
!		2003 FORMAT (15F10.5)
			
	DEALLOCATE(Autocor_helper)
	!Performing Crosscorrelation
	sumX = 0
	sumX2 = 0
	sumY = 0 
	sumY2 = 0 
	sumXY = 0
	offset = 0
	helper = 0


	ALLOCATE (Crosscor_helper(segment_numbers, Correlation_steps))
	Crosscor_helper = 0

	outer_cross : DO i=1, segment_numbers
		helper = segment_rows * (i-1)
		middle_cross : DO j=1, Correlation_steps
			inner_cross : DO k = 1, boundary
				sumX = sumX + Data_array_cut(2 , (k + 199 + helper))
				sumX2 = sumX2 + Data_array_cut(2 , (k +199 + helper))* &
						Data_array_cut(2 , (k +199 + helper)) 
				sumY = sumY + Data_array_cut(3, (k+ helper+offset))
				sumY2 = sumY2 + Data_array_cut(3, (k+ helper+offset))* &
						Data_array_cut(3, (k+ helper+offset))
				sumXY = sumXY + Data_array_cut(2, (k + 199 + helper)) * &
						Data_array_cut(3, (k+ helper+offset))

			END DO inner_cross
		
			Crosscor_helper(i,j) = (boundary*sumXY - sumX * sumY)/(SQRT(boundary*sumX2- &
									sumX*sumX)* SQRT(boundary*sumY2-sumY*sumY))
			! Adjust variables
			offset=offset +1
			sumX = 0
			sumX2 = 0
			sumY = 0 
			sumY2 = 0 
			sumXY = 0 			

		END DO middle_cross

		offset = 0		!zero offset helper
	END DO outer_cross

	!average the autocorrelation-segments and add time
	DO i = 1, Correlation_steps
			Crosscorrelation_raw (1, i) = -REAL(199)/frequency + real((i-1))/frequency
			
			DO j = 1, segment_numbers
				average_value = average_value + Crosscor_helper (j , i)
			END DO
					
			Crosscorrelation_raw (2, i) = average_value/segment_numbers
			average_value = 0
	END DO
	 
		WRITE (*,*) 'Finished calculation of crosscorrelation_raw! '

!		OPEN (UNIT = 2000, FILE = 'cross.txt', STATUS = 'REPLACE', ACTION ='WRITE', IOSTAT=status)
!		WRITE (2000, 2001) ((Crosscor_helper (i,j), i=1,segment_numbers), j=1, Correlation_steps)
!		2001 FORMAT (15F10.5)
			
	DEALLOCATE(Crosscor_helper)
	
		
END SUBROUTINE Correlation_raw


!_______________________________________________________________________________________________________________


SUBROUTINE Correlation_scales_raw (Autocorrelation_raw, Crosscorrelation_raw, Correlation_steps, &
		   CorScale_results_raw, frequency, delta_x, velocity)
IMPLICIT NONE

	INTEGER, PARAMETER :: DBL = SELECTED_REAL_KIND(p=15)	!Precision of the Real variables
	INTEGER, INTENT (IN) :: Correlation_steps
	INTEGER, INTENT (IN) :: frequency
	REAL, INTENT (IN) :: delta_x
	REAL(kind=DBL), INTENT (INOUT) :: velocity
	REAL(kind=DBL), INTENT (IN) :: Autocorrelation_raw (2, Correlation_steps)
	REAL(kind=DBL), INTENT (IN) :: Crosscorrelation_raw (2, Correlation_steps)
	REAL(kind=DBL), INTENT (OUT) ::CorScale_results_raw(15)
	INTEGER :: i = 0 , j = 0, counter, locator
	REAL(kind=DBL) :: helper
	INTEGER :: helper_int
	REAL(kind=DBL) :: min_auto, min_auto_loc, min_cross, min_cross_loc 
	REAL(kind=DBL) :: T_point_five
	REAL(kind=DBL) :: Zero_crossing_auto
	REAL(kind=DBL) :: Txx
	REAL(kind=DBL) :: Rxz_max
	REAL(kind=DBL) :: T_Rxz_max
	REAL(kind=DBL) :: T_half_Rxz_max 
	REAL(kind=DBL) :: Zero_crossing_cross 
	REAL(kind=DBL) :: Txz
	REAL(kind=DBL) :: Tu 
	REAL(kind=DBL) :: v_prime

	velocity = 0
	CorScale_results_raw = 0
	helper = 0
	T_point_five = 0
	Zero_crossing_auto =0
	min_auto_loc = 0
	min_auto = 0
	Txx = 0
	Rxz_max = 0
	T_Rxz_max = 0
	Zero_crossing_cross = 0
	min_cross_loc = 0
	min_cross = 0
	T_half_Rxz_max = 0
	Txz = 0
	Tu = 0
	Velocity = 0
	counter=0
	locator = 0
	helper_int = 0
	v_prime = 0

	!Autocorrelation_results
	!Calculation of Time when Autocorrelation is 0.5
	DO i=1, Correlation_steps
		helper = Autocorrelation_raw(2, i)
		IF (helper == 0.5) THEN
			T_point_five = Autocorrelation_raw(1,i)
		EXIT
			ELSE IF (helper < 0.5) THEN
			    T_point_five = Autocorrelation_raw(1,i-1) + (Autocorrelation_raw(1,i) - &
				Autocorrelation_raw(1,i-1)) / (Autocorrelation_raw(2,i) - Autocorrelation_raw(2,i-1))&
				 * (0.5-Autocorrelation_raw(2,i-1))	
			EXIT
		END IF
	END DO

	helper=0
	counter = 0

	!Calculation of Zero crossing of Autocorrelation function
	DO i=1, Correlation_steps
		helper = Autocorrelation_raw(2,i)
		Zero_crossing_auto = 0
		counter = counter +1
		IF (helper == 0) THEN
			Zero_crossing_auto = Autocorrelation_raw(1,i)
				EXIT
			ELSE IF (helper < 0.) THEN
				Zero_crossing_auto = Autocorrelation_raw(1,i-1) + (Autocorrelation_raw(1,i) - &
				Autocorrelation_raw(1,i-1)) / (Autocorrelation_raw(2,i) - Autocorrelation_raw(2,i-1)) &
				* (0-Autocorrelation_raw(2,i-1))
				EXIT
				ELSE 
				Zero_crossing_auto = 99	
		END IF
	END DO	

	!Calculation of Txx
	min_auto=1
	locator = 0
	min_auto_loc = 0
	Txx =0
	IF (INT(Zero_crossing_auto) == 99) THEN
		DO j = 1, Correlation_steps-1
			min_auto = Min(Autocorrelation_raw(2,j), min_auto)
		END DO

		DO i = 1, Correlation_steps-1
		locator = locator +1
			IF (Autocorrelation_raw(2,i) - min_auto < 0.00001) THEN
			min_auto_loc = Autocorrelation_raw(1,locator)
			EXIT
			END IF
		END DO

		DO i = 1, locator
			IF (i ==1000) THEN
			Txx = Txx + (Autocorrelation_raw(2,i)/frequency)
			ELSE
			Txx = Txx + (Autocorrelation_raw(2,i)+Autocorrelation_raw(2,i+1))/(2*frequency)
			END IF
		END DO
	ELSE 
		DO j = 1, counter
			IF (j ==1000) THEN
			Txx = Txx + (Autocorrelation_raw(2,j)/frequency)
			ELSE
			Txx = Txx + (Autocorrelation_raw(2,j)+Autocorrelation_raw(2,j+1))/(2*frequency)
			END IF
		END DO 
		min_auto=0
		min_auto_loc=Zero_crossing_auto
	END IF
		
	helper =0
	counter =0
	Rxz_max = 0
	!Crosscorrelation_results
	!Find maximum Crosscorrelation value (Rxz)max
	DO i = 1, Correlation_steps
		helper = Crosscorrelation_raw(2,i)
		Rxz_max = Max(Rxz_max,helper)
	END DO	
		
	helper = 0 
	locator = 0
	T_Rxz_max = 0

	!Corresponding Time for Rxz_max
	DO i = 1, Correlation_steps
		locator = locator +1
		IF (Rxz_max - Crosscorrelation_raw(2,i) < 0.0000001) THEN
		T_Rxz_max = Crosscorrelation_raw(1,locator)
		EXIT 
		END IF
	END DO

	T_half_Rxz_max = 0

	!Calculation of Time for 0.5*Rxz_max
	DO i=locator, Correlation_steps
		IF (Crosscorrelation_raw (2,i) <= 0.5*Rxz_max) THEN
			T_half_Rxz_max = Crosscorrelation_raw(1,i-1) + (Crosscorrelation_raw(1,i) - &
							Crosscorrelation_raw(1,i-1)) / (Crosscorrelation_raw(2,i) - &
							Crosscorrelation_raw(2,i-1)) * (0.5*Rxz_max-Crosscorrelation_raw(2,i-1))
			EXIT
		END IF
	END DO

	counter =-1
	helper = 0

	!Time where the values of Rxz_max cross the x-axis for the first time right of Rxy_max
	DO i=locator, Correlation_steps
		helper = Crosscorrelation_raw(2,i)
		Zero_crossing_cross = 0
		counter = counter +1
		IF (helper == 0) THEN
			Zero_crossing_cross = Crosscorrelation_raw(1,i)
				EXIT
			ELSE IF (helper < 0.) THEN
				Zero_crossing_cross = Crosscorrelation_raw(1,i-1) + (Crosscorrelation_raw(1,i) - &
			    Crosscorrelation_raw(1,i-1))  / (Crosscorrelation_raw(2,i) - Crosscorrelation_raw(2,i-1))&
				* (0-Crosscorrelation_raw(2,i-1))
				EXIT
				ELSE 
				Zero_crossing_cross = 98	
		END IF
	END DO	

	!Calculation of Txz
	min_cross = Rxz_max
	min_cross_loc = 0
	Txz = 0
	IF (INT(Zero_crossing_cross) == 98) THEN
		DO j = locator, Correlation_steps-1
			min_cross = Min(Crosscorrelation_raw(2,j), min_cross)
		END DO
		
		helper_int = 0
		DO i = locator, Correlation_steps-1
		helper_int = helper_int +1
			IF (Crosscorrelation_raw(2,i) - min_cross < 0.0000001) THEN
			min_cross_loc = Crosscorrelation_raw(1,(helper_int+locator-1))
			EXIT
			END IF
		END DO

		DO i = locator, (helper_int +locator-1)
			IF (i ==1000) THEN
			Txz = Txz + (Crosscorrelation_raw(2,i)/frequency)
			ELSE
			Txz = Txz + (Crosscorrelation_raw(2,i)+Crosscorrelation_raw(2,i+1))/(2*frequency)
			END IF
		END DO
	ELSE 
		DO j = locator, (locator+counter)
			IF (j ==1000) THEN
			Txz = Txz + (Crosscorrelation_raw(2,j)/frequency)
			ELSE
			Txz = Txz + (Crosscorrelation_raw(2,j)+Crosscorrelation_raw(2,j+1))/(2*frequency)
			END IF
		END DO 
		min_cross = 0
		min_cross_loc = Zero_crossing_cross
	END IF
		
	!Calculation of turbulence intensity
	IF (((T_half_Rxz_max - T_Rxz_max)*(T_half_Rxz_max - T_Rxz_max)- &
		 T_point_five*T_point_five) <0) THEN
		 Tu = 97
	ELSE
	Tu = 0.851 * SQRT((T_half_Rxz_max - T_Rxz_max)*(T_half_Rxz_max - T_Rxz_max)- &
		 T_point_five*T_point_five)/T_Rxz_max
	END IF
		
	!Calculation of time averaged interfacial velocity
	velocity = delta_x/1000/T_Rxz_max	

	v_prime = Tu * velocity

	!store results into array CorScale_results  :T_point_five, Zero_crossing_auto, min_auto_loc, min_auto,Txx, &
						!T_Rxz_max, T_half_Rxz_max, Zero_crossing_cross, min_cross_loc, min_cross,Txz, Tu, Velocity 
	CorScale_results_raw = (/ T_point_five, Zero_crossing_auto, min_auto_loc, min_auto, Txx, Rxz_max, T_Rxz_max, &
						Zero_crossing_cross, min_cross_loc, min_cross, T_half_Rxz_max, Txz, Tu, Velocity, v_prime /)


END SUBROUTINE Correlation_scales_raw

!_______________________________________________________________________________________________________________

!______________________________________________________________________________________________________

SUBROUTINE Correlation_high (Filtered_array_high_cut, Filtered_array_mid_cut, number_columns, number_rows_cut, Correlation_steps, &  
			                Autocorrelation_high, Crosscorrelation_high, frequency, segment_numbers, Auto_fac, Cross_fac)
IMPLICIT NONE

	INTEGER, PARAMETER :: DBL = SELECTED_REAL_KIND(p=15)	!Precision of the Real variables
	INTEGER, INTENT (IN) :: number_columns
	INTEGER, INTENT (IN) :: number_rows_cut
	INTEGER, INTENT (IN) :: Correlation_steps
	INTEGER, INTENT (IN) :: frequency
	INTEGER, INTENT (IN) :: segment_numbers
	REAL(kind=DBL), INTENT (IN) :: Filtered_array_high_cut (number_columns, number_rows_cut)
	REAL(kind=DBL), INTENT (IN) :: Filtered_array_mid_cut (number_columns, number_rows_cut)
	REAL(kind=DBL), INTENT (OUT) :: Autocorrelation_high (2, Correlation_steps)
	REAL(kind=DBL), INTENT (OUT) :: Crosscorrelation_high(2, Correlation_steps)
	REAL(kind=DBL), INTENT(OUT) :: Auto_fac
	REAL(kind=DBL), INTENT(OUT) :: Cross_fac
	INTEGER :: segment_rows							!number of rows in each segment
	INTEGER :: helper								!helper for faster correlations
	INTEGER :: boundary								!helper for faster correlation
	INTEGER :: i =0, j = 0, k=0							!integer for loop
	REAL(kind=DBL), ALLOCATABLE, DIMENSION  (:,:) :: Autocor_helper	
	REAL(kind=DBL), ALLOCATABLE, DIMENSION  (:,:) :: Crosscor_helper	
	REAL(kind=DBL) :: sumX, sumY, sumX2, sumY2, sumXY !Variables for correl
	INTEGER :: offset												 !integer for offsetting
	INTEGER :: status
	REAL(kind=DBL) :: average_value				!Helper to calculate the average values 
	
	! calculation for factor
	REAL(kind=DBL) :: Crosscor_factor(segment_numbers)
	REAL(kind=DBL) :: Autocor_factor(segment_numbers)
	REAL(kind=DBL) :: Boundary_factor
	REAL(kind=DBL) :: average_factor
	REAL(kind=DBL), ALLOCATABLE, DIMENSION  (:,:) :: denominator_array	
!	REAL(kind=DBL), ALLOCATABLE, DIMENSION  (:,:) :: Crosscor_factor_helper	
!	REAL(kind=DBL), ALLOCATABLE, DIMENSION  (:,:) :: Autocor_factor_helper	
													
	
	Autocorrelation_high = 0
	Crosscorrelation_high = 0
	segment_rows	=0
	helper = 0
	boundary = 0
	sumX = 0
	sumY = 0
	sumX2 = 0
	sumY2 = 0
	sumXY = 0
	offset = 0
	status = 0
	average_value =0
	Crosscor_factor = 0 
	Autocor_factor = 0 
	Boundary_factor = 0
	average_factor = 0

	segment_rows = number_rows_cut/segment_numbers
	boundary = segment_rows - Correlation_steps	

	Boundary_factor = number_rows_cut - Correlation_steps

	!FACTOR calculations

	ALLOCATE (denominator_array(number_columns, number_rows_cut))

		Do j = 1, number_rows_cut
			denominator_array (2,j) =  (Filtered_array_mid_cut (2, j) + Filtered_array_high_cut (2, j)) * (Filtered_array_mid_cut (2, j) + Filtered_array_high_cut (2, j))
			denominator_array (3,j) =  (Filtered_array_mid_cut (3, j) + Filtered_array_high_cut (3, j)) * (Filtered_array_mid_cut (3, j) + Filtered_array_high_cut (3, j))
		END DO

	! cross-correlation -> factor


	DO i=1, segment_numbers
		helper = segment_rows * (i-1)
			DO k = 1, boundary
				sumX = sumX + (Filtered_array_high_cut(2 , (k + 199+ helper))) * (Filtered_array_high_cut(2 , (k + 199+ helper)))
				sumY = sumY + (Filtered_array_high_cut(3, (k+ helper))) * (Filtered_array_high_cut(3, (k+ helper)))
				
				sumX2 = sumX2 + denominator_array(2, (k +199+ helper))
				sumY2 = sumY2 + denominator_array(3, (k + helper))

			END DO
			Crosscor_factor (i) = (SQRT(sumX * sumY)) / (SQRT(sumX2 * sumY2))
			sumX = 0
			sumX2 = 0
			sumY = 0 
			sumY2 = 0 
	END DO

			! Adjust variables



!	ALLOCATE (Crosscor_factor_helper(segment_numbers, Correlation_steps))
!	Crosscor_factor_helper = 0

!	outer_cross : DO i=1, segment_numbers
!		helper = segment_rows * (i-1)
!		middle_cross : DO j=1, Correlation_steps
!			inner_cross : DO k = 1, boundary
!				sumX = sumX + (Filtered_array_high_cut(2 , (k + 199 + helper))) * (Filtered_array_high_cut(2 , (k + 199 + helper)))
!				sumY = sumY + (Filtered_array_high_cut(3, (k+ helper+offset))) * (Filtered_array_high_cut(3, (k+ helper+offset)))
				
!				sumX2 = sumX2 + denominator_array(2, (k +199 + helper))
!				sumY2 = sumY2 + denominator_array(3, (k + helper + offset))

!			END DO inner_cross
!			Crosscor_factor_helper(i,j) = (SQRT(sumX * sumY)) / (SQRT(sumX2 * sumY2))

			! Adjust variables
!			offset=offset +1
!			sumX = 0
!			sumX2 = 0
!			sumY = 0 
!			sumY2 = 0 

!		END DO middle_cross

!		offset = 0		!zero offset helper
!	END DO outer_cross

!	sumX = 0
!	sumX2 = 0
!	sumY = 0 
!	sumY2 = 0 
!	offset = 0
!	helper = 0


	! Auto-correlation -> factor
	DO i=1, segment_numbers
		helper = segment_rows * (i-1)
			DO k = 1, boundary
				sumX = sumX + (Filtered_array_high_cut(2 , (k + helper))) * (Filtered_array_high_cut(2 , (k + helper)))
	!			sumY = sumY + (Filtered_array_high_cut(2, (k + helper))) * (Filtered_array_high_cut(2, (k + helper)))
				
				sumX2 = sumX2 + denominator_array(2, (k + helper))
	!			sumY2 = sumY2 + denominator_array(2, (k + helper))

			END DO
			Autocor_factor(i) = sumX / sumX2
			sumX = 0
			sumX2 = 0
			sumY = 0 
			sumY2 = 0 
	END DO
			


!	ALLOCATE (Autocor_factor_helper(segment_numbers, Correlation_steps))
!	Autocor_factor_helper = 0

!	DO i=1, segment_numbers
!		helper = segment_rows * (i-1)
!		DO j=1, Correlation_steps
!			DO k = 1, boundary
!				sumX = sumX + (Filtered_array_high_cut(2 , (k + helper))) * (Filtered_array_high_cut(2 , (k + helper)))
!				sumY = sumY + (Filtered_array_high_cut(2, (k+ helper+offset))) * (Filtered_array_high_cut(2, (k+ helper+offset)))
				
!				sumX2 = sumX2 + denominator_array(2, (k + helper))
!				sumY2 = sumY2 + denominator_array(2, (k + helper + offset))

!			END DO 
!			Autocor_factor_helper(i,j) = (SQRT(sumX * sumY)) / (SQRT(sumX2 * sumY2))

			! Adjust variables
!			offset=offset +1
!			sumX = 0
!			sumX2 = 0
!			sumY = 0 
!			sumY2 = 0 

!		END DO 

!		offset = 0		!zero offset helper
!	END DO 


!	sumX = 0
!	sumX2 = 0
!	sumY = 0 
!	sumY2 = 0 
!	offset = 0
!	helper = 0
!

	!Performing Auto-correlation
	ALLOCATE (Autocor_helper(segment_numbers, Correlation_steps))
	Autocor_helper = 0			
	outer : DO i=1, segment_numbers
			helper = segment_rows * (i-1)
		middle : DO j=1, Correlation_steps
			inner : DO k = 1, boundary
					sumX = sumX + Filtered_array_high_cut(2 , (k + helper))
					sumX2 = sumX2 + Filtered_array_high_cut(2 , (k+ helper))* &
							Filtered_array_high_cut(2 , (k+ helper)) 
					sumY = sumY + Filtered_array_high_cut(2, (k+ helper+offset))
					sumY2 = sumY2 + Filtered_array_high_cut(2, (k+ helper+offset)) &
							*Filtered_array_high_cut(2, (k+ helper+offset))
					sumXY = sumXY + Filtered_array_high_cut(2,(k+ helper)) * &
							Filtered_array_high_cut(2,(k+ helper+offset))
			END DO inner
		
			Autocor_helper(i,j) = (sumXY - sumX * sumY/boundary) /(SQRT(sumX2-sumX*sumX/ &
								  boundary) * SQRT(sumY2-sumY*sumY/boundary))
			! Adjust variables
			offset=offset +1
			sumX = 0
			sumX2 = 0
			sumY = 0 
			sumY2 = 0 
			sumXY = 0 			

		END DO middle
		offset = 0		!zero offset helper
	END DO outer
	helper = 0

	!average the autocorrelation-segments and add time
	DO i = 1, Correlation_steps
			Autocorrelation_high (1, i) = real((i-1))/frequency
			
			DO j = 1, segment_numbers
					average_value = average_value + (Autocor_helper (j , i) * Autocor_factor(j))
					
			END DO 
				
				Autocorrelation_high (2, i) = average_value/segment_numbers
				average_value = 0
				 
				
	END DO

	DO j = 1, segment_numbers
	average_factor = average_factor + Autocor_factor(j)
	END Do
	Auto_fac = average_factor/segment_numbers

	average_factor = 0

	WRITE (*,*) 'Finished calculation of autocorrelation_high! '

!		OPEN (UNIT = 2002, FILE = 'auto.txt', STATUS = 'REPLACE', ACTION ='WRITE', IOSTAT=status)
!		WRITE (2002, 2003) ((Autocor_helper (i,j), i=1,segment_numbers), j=1, Correlation_steps)
!		2003 FORMAT (15F10.5)
			
	DEALLOCATE(Autocor_helper)
!	DEALLOCATE(Autocor_factor_helper)
	!Performing Crosscorrelation
	sumX = 0
	sumX2 = 0
	sumY = 0 
	sumY2 = 0 
	sumXY = 0
	offset = 0
	helper = 0


	ALLOCATE (Crosscor_helper(segment_numbers, Correlation_steps))
	Crosscor_helper = 0

	DO i=1, segment_numbers
		helper = segment_rows * (i-1)
		DO j=1, Correlation_steps
			DO k = 1, boundary
				sumX = sumX + Filtered_array_high_cut(2 , (k + 199 + helper))
				sumX2 = sumX2 + Filtered_array_high_cut(2 , (k +199 + helper))* &
						Filtered_array_high_cut(2 , (k +199 + helper)) 
				sumY = sumY + Filtered_array_high_cut(3, (k+ helper+offset))
				sumY2 = sumY2 + Filtered_array_high_cut(3, (k+ helper+offset))* &
						Filtered_array_high_cut(3, (k+ helper+offset))
				sumXY = sumXY + Filtered_array_high_cut(2, (k + 199 + helper)) * &
						Filtered_array_high_cut(3, (k+ helper+offset))

			END DO 
		
			Crosscor_helper(i,j) = (boundary*sumXY - sumX * sumY)/(SQRT(boundary*sumX2- &
									sumX*sumX)* SQRT(boundary*sumY2-sumY*sumY))
			! Adjust variables
			offset=offset +1
			sumX = 0
			sumX2 = 0
			sumY = 0 
			sumY2 = 0 
			sumXY = 0 			

		END DO

		offset = 0		!zero offset helper
	END DO

	!average the autocorrelation-segments and add time
	DO i = 1, Correlation_steps
			Crosscorrelation_high (1, i) = -REAL(199)/frequency + real((i-1))/frequency
			
			DO j = 1, segment_numbers
				average_value = average_value + (Crosscor_helper (j , i) * Crosscor_factor(j))
			END DO
					
			Crosscorrelation_high (2, i) = average_value/segment_numbers
			average_value = 0
	END DO

	DO j = 1, segment_numbers
	average_factor = average_factor + Crosscor_factor(j)
	END Do
	Cross_fac = average_factor/segment_numbers

	average_factor = 0
	 
		WRITE (*,*) 'Finished calculation of crosscorrelation_high! '

!		OPEN (UNIT = 2000, FILE = 'cross.txt', STATUS = 'REPLACE', ACTION ='WRITE', IOSTAT=status)
!		WRITE (2000, 2001) ((Crosscor_helper (i,j), i=1,segment_numbers), j=1, Correlation_steps)
!		2001 FORMAT (15F10.5)
			
	DEALLOCATE(Crosscor_helper)
!	DEALLOCATE(Crosscor_factor_helper)
	DEALLOCATE(denominator_array)
	
		
END SUBROUTINE Correlation_high


!_______________________________________________________________________________________________________________


SUBROUTINE Correlation_scales_high (Autocorrelation_high, Crosscorrelation_high, Correlation_steps, &
		   CorScale_results_high, frequency, delta_x, velocity, Auto_fac, Cross_fac)
IMPLICIT NONE

	INTEGER, PARAMETER :: DBL = SELECTED_REAL_KIND(p=15)	!Precision of the Real variables
	INTEGER, INTENT (IN) :: Correlation_steps
	INTEGER, INTENT (IN) :: frequency
	REAL, INTENT (IN) :: delta_x
	REAL(kind=DBL), INTENT (INOUT) :: velocity
	REAL(kind=DBL), INTENT (INOUT) :: Auto_fac
	REAL(kind=DBL), INTENT (INOUT) :: Cross_fac
	REAL(kind=DBL), INTENT (IN) :: Autocorrelation_high (2, Correlation_steps)
	REAL(kind=DBL), INTENT (IN) :: Crosscorrelation_high (2, Correlation_steps)
	REAL(kind=DBL), INTENT (OUT) ::CorScale_results_high(18)
	INTEGER :: i = 0 , j = 0, counter, locator
	REAL(kind=DBL) :: helper
	INTEGER :: helper_int
	REAL(kind=DBL) :: min_auto, min_auto_loc, min_cross, min_cross_loc 
	REAL(kind=DBL) :: T_point_five
	REAL(kind=DBL) :: Zero_crossing_auto
	REAL(kind=DBL) :: Txx
	REAL(kind=DBL) :: Rxz_max
	REAL(kind=DBL) :: T_Rxz_max
	REAL(kind=DBL) :: T_half_Rxz_max 
	REAL(kind=DBL) :: Zero_crossing_cross 
	REAL(kind=DBL) :: Txz
	REAL(kind=DBL) :: Tu 
	REAL(kind=DBL) :: v_prime
	REAL(kind=DBL) :: Rxx_max
	REAL(kind=DBL) :: helper_RxxMax_half

	velocity = 0
	CorScale_results_high = 0
	helper = 0
	T_point_five = 0
	Zero_crossing_auto =0
	min_auto_loc = 0
	min_auto = 0
	Txx = 0
	Rxz_max = 0
	T_Rxz_max = 0
	Zero_crossing_cross = 0
	min_cross_loc = 0
	min_cross = 0
	T_half_Rxz_max = 0
	Txz = 0
	Tu = 0
	Velocity = 0
	counter=0
	locator = 0
	helper_int = 0
	v_prime = 0
	Rxx_max = 0

	
	Rxx_max  = Autocorrelation_high (2,1)
	helper_RxxMax_half = Rxx_max/2

	!Autocorrelation_results
	!Calculation of Time when Autocorrelation is 0.5*Rxx_max
	DO i=1, Correlation_steps
		helper = Autocorrelation_high(2, i)
		IF (helper == helper_RxxMax_half) THEN
			T_point_five = Autocorrelation_high(1,i)
		EXIT
			ELSE IF (helper < helper_RxxMax_half) THEN
			    T_point_five = Autocorrelation_high(1,i-1) + (Autocorrelation_high(1,i) - &
				Autocorrelation_high(1,i-1)) / (Autocorrelation_high(2,i) - Autocorrelation_high(2,i-1))&
				 * (helper_RxxMax_half-Autocorrelation_high(2,i-1))	
			EXIT
		END IF
	END DO

	helper=0
	counter = 0

	!Calculation of Zero crossing of Autocorrelation function
	DO i=1, Correlation_steps
		helper = Autocorrelation_high(2,i)
		Zero_crossing_auto = 0
		counter = counter +1
		IF (helper == 0) THEN
			Zero_crossing_auto = Autocorrelation_high(1,i)
				EXIT
			ELSE IF (helper < 0.) THEN
				Zero_crossing_auto = Autocorrelation_high(1,i-1) + (Autocorrelation_high(1,i) - &
				Autocorrelation_high(1,i-1)) / (Autocorrelation_high(2,i) - Autocorrelation_high(2,i-1)) &
				* (0-Autocorrelation_high(2,i-1))
				EXIT
				ELSE 
				Zero_crossing_auto = 99	
		END IF
	END DO	

	!Calculation of Txx
	min_auto=1
	locator = 0
	min_auto_loc = 0
	Txx =0
	IF (INT(Zero_crossing_auto) == 99) THEN
		DO j = 1, Correlation_steps-1
			min_auto = Min(Autocorrelation_high(2,j), min_auto)
		END DO

		DO i = 1, Correlation_steps-1
		locator = locator +1
			IF (Autocorrelation_high(2,i) - min_auto < 0.00001) THEN
			min_auto_loc = Autocorrelation_high(1,locator)
			EXIT
			END IF
		END DO

		DO i = 1, locator
			IF (i ==1000) THEN
			Txx = Txx + (Autocorrelation_high(2,i)/frequency)
			ELSE
			Txx = Txx + (Autocorrelation_high(2,i)+Autocorrelation_high(2,i+1))/(2*frequency)
			END IF
		END DO
	ELSE 
		DO j = 1, counter
			IF (j ==1000) THEN
			Txx = Txx + (Autocorrelation_high(2,j)/frequency)
			ELSE
			Txx = Txx + (Autocorrelation_high(2,j)+Autocorrelation_high(2,j+1))/(2*frequency)
			END IF
		END DO 
		min_auto=0
		min_auto_loc=Zero_crossing_auto
	END IF
		
	helper =0
	counter =0
	Rxz_max = 0
	!Crosscorrelation_results
	!Find maximum Crosscorrelation value (Rxz)max
	DO i = 1, Correlation_steps
		helper = Crosscorrelation_high(2,i)
		Rxz_max = Max(Rxz_max,helper)
	END DO	
		
	helper = 0 
	locator = 0
	T_Rxz_max = 0

	!Corresponding Time for Rxz_max
	DO i = 1, Correlation_steps
		locator = locator +1
		IF (Rxz_max - Crosscorrelation_high(2,i) < 0.0000001) THEN
		T_Rxz_max = Crosscorrelation_high(1,locator)
		EXIT 
		END IF
	END DO

	T_half_Rxz_max = 0

	!Calculation of Time for 0.5*Rxz_max
	DO i=locator, Correlation_steps
		IF (Crosscorrelation_high (2,i) <= 0.5*Rxz_max) THEN
			T_half_Rxz_max = Crosscorrelation_high(1,i-1) + (Crosscorrelation_high(1,i) - &
							Crosscorrelation_high(1,i-1)) / (Crosscorrelation_high(2,i) - &
							Crosscorrelation_high(2,i-1)) * (0.5*Rxz_max-Crosscorrelation_high(2,i-1))
			EXIT
		END IF
	END DO

	counter =-1
	helper = 0

	!Time where the values of Rxz_max cross the x-axis for the first time right of Rxy_max
	DO i=locator, Correlation_steps
		helper = Crosscorrelation_high(2,i)
		Zero_crossing_cross = 0
		counter = counter +1
		IF (helper == 0) THEN
			Zero_crossing_cross = Crosscorrelation_high(1,i)
				EXIT
			ELSE IF (helper < 0.) THEN
				Zero_crossing_cross = Crosscorrelation_high(1,i-1) + (Crosscorrelation_high(1,i) - &
			    Crosscorrelation_high(1,i-1))  / (Crosscorrelation_high(2,i) - Crosscorrelation_high(2,i-1))&
				* (0-Crosscorrelation_high(2,i-1))
				EXIT
				ELSE 
				Zero_crossing_cross = 98	
		END IF
	END DO	

	!Calculation of Txz
	min_cross = Rxz_max
	min_cross_loc = 0
	Txz = 0
	IF (INT(Zero_crossing_cross) == 98) THEN
		DO j = locator, Correlation_steps-1
			min_cross = Min(Crosscorrelation_high(2,j), min_cross)
		END DO
		
		helper_int = 0
		DO i = locator, Correlation_steps-1
		helper_int = helper_int +1
			IF (Crosscorrelation_high(2,i) - min_cross < 0.0000001) THEN
			min_cross_loc = Crosscorrelation_high(1,(helper_int+locator-1))
			EXIT
			END IF
		END DO

		DO i = locator, (helper_int +locator-1)
			IF (i ==1000) THEN
			Txz = Txz + (Crosscorrelation_high(2,i)/frequency)
			ELSE
			Txz = Txz + (Crosscorrelation_high(2,i)+Crosscorrelation_high(2,i+1))/(2*frequency)
			END IF
		END DO
	ELSE 
		DO j = locator, (locator+counter)
			IF (j ==1000) THEN
			Txz = Txz + (Crosscorrelation_high(2,j)/frequency)
			ELSE
			Txz = Txz + (Crosscorrelation_high(2,j)+Crosscorrelation_high(2,j+1))/(2*frequency)
			END IF
		END DO 
		min_cross = 0
		min_cross_loc = Zero_crossing_cross
	END IF
		
	!Calculation of turbulence intensity
	IF (((T_half_Rxz_max - T_Rxz_max)*(T_half_Rxz_max - T_Rxz_max)- &
		 T_point_five*T_point_five) <0) THEN
		 Tu = 97
	ELSE
	Tu = 0.851 * SQRT((T_half_Rxz_max - T_Rxz_max)*(T_half_Rxz_max - T_Rxz_max)- &
		 T_point_five*T_point_five)/T_Rxz_max
	END IF
		
	!Calculation of time averaged interfacial velocity
	velocity = delta_x/1000/T_Rxz_max
	
	v_prime = Tu * velocity

	!store results into array CorScale_results  :T_point_five, Zero_crossing_auto, min_auto_loc, min_auto,Txx, &
						!T_Rxz_max, T_half_Rxz_max, Zero_crossing_cross, min_cross_loc, min_cross,Txz, Tu, Velocity, Rxx_max 
	CorScale_results_high = (/ T_point_five, Zero_crossing_auto, min_auto_loc, min_auto, Txx, Rxz_max, T_Rxz_max, &
						Zero_crossing_cross, min_cross_loc, min_cross, T_half_Rxz_max, Txz, Tu, Velocity, v_prime, Rxx_max, Auto_fac, Cross_fac /)


END SUBROUTINE Correlation_scales_high

!____________________________________
!______________________________________________________________________________________________________

SUBROUTINE Correlation_mid (Filtered_array_mid_cut, Filtered_array_high_cut, number_columns, number_rows_cut, Correlation_steps, &  
			                Autocorrelation_mid, Crosscorrelation_mid, frequency, segment_numbers, Auto_fac, Cross_fac)
IMPLICIT NONE

	INTEGER, PARAMETER :: DBL = SELECTED_REAL_KIND(p=15)	!Precision of the Real variables
	INTEGER, INTENT (IN) :: number_columns
	INTEGER, INTENT (IN) :: number_rows_cut
	INTEGER, INTENT (IN) :: Correlation_steps
	INTEGER, INTENT (IN) :: frequency
	INTEGER, INTENT (IN) :: segment_numbers
	REAL(kind=DBL), INTENT(OUT) :: Auto_fac
	REAL(kind=DBL), INTENT(OUT) :: Cross_fac
	REAL(kind=DBL), INTENT (IN) :: Filtered_array_mid_cut (number_columns, number_rows_cut)
	REAL(kind=DBL), INTENT (IN) :: Filtered_array_high_cut (number_columns, number_rows_cut)
	REAL(kind=DBL), INTENT (OUT) :: Autocorrelation_mid (2, Correlation_steps)
	REAL(kind=DBL), INTENT (OUT) :: Crosscorrelation_mid(2, Correlation_steps)
	INTEGER :: segment_rows							!number of rows in each segment
	INTEGER :: helper								!helper for faster correlations
	INTEGER :: boundary								!helper for faster correlation
	INTEGER :: i =0, j = 0, k=0							!integer for loop
	REAL(kind=DBL), ALLOCATABLE, DIMENSION  (:,:) :: Autocor_helper	
	REAL(kind=DBL), ALLOCATABLE, DIMENSION  (:,:) :: Crosscor_helper	
	REAL(kind=DBL) :: sumX, sumY, sumX2, sumY2, sumXY !Variables for correl
	INTEGER :: offset												 !integer for offsetting
	INTEGER :: status
	REAL(kind=DBL) :: average_value				!Helper to calculate the average values 
	REAL(kind=DBL) :: average_factor	
	
	! calculation for factor
	REAL(kind=DBL), ALLOCATABLE, DIMENSION  (:,:) :: denominator_array	
	REAL(kind=DBL), ALLOCATABLE, DIMENSION  (:,:) :: Crosscor_factor_helper	
	REAL(kind=DBL), ALLOCATABLE, DIMENSION  (:,:) :: Autocor_factor_helper	

	REAL(kind=DBL) :: Crosscor_factor(segment_numbers)
	REAL(kind=DBL) :: Autocor_factor(segment_numbers)
	REAL(kind=DBL) :: Boundary_factor
													
	
	Autocorrelation_mid = 0
	Crosscorrelation_mid = 0
	segment_rows	=0
	helper = 0
	boundary = 0
	sumX = 0
	sumY = 0
	sumX2 = 0
	sumY2 = 0
	sumXY = 0
	offset = 0
	status = 0
	average_value =0
	Crosscor_factor = 0 
	Autocor_factor = 0 
	Boundary_factor = 0
	average_factor = 0

	segment_rows = number_rows_cut/segment_numbers
	boundary = segment_rows - Correlation_steps	

	Boundary_factor = number_rows_cut - Correlation_steps	

	!FACTOR calculations
	ALLOCATE (denominator_array(number_columns, number_rows_cut))
	denominator_array = 0

		Do j = 1, number_rows_cut
			denominator_array (2,j) =  (Filtered_array_mid_cut (2, j) + Filtered_array_high_cut (2, j)) * (Filtered_array_mid_cut (2, j) + Filtered_array_high_cut (2, j))
			denominator_array (3,j) =  (Filtered_array_mid_cut (3, j) + Filtered_array_high_cut (3, j)) * (Filtered_array_mid_cut (3, j) + Filtered_array_high_cut (3, j))
		END DO

	! cross-correlation -> factor
	DO i=1, segment_numbers
		helper = segment_rows * (i-1)
			DO k = 1, boundary
				sumX = sumX + (Filtered_array_mid_cut(2 , (k + 199 + helper))) * (Filtered_array_mid_cut(2 , (k + 199 + helper)))
				sumY = sumY + (Filtered_array_mid_cut(3, (k + helper))) * (Filtered_array_mid_cut(3, (k + helper)))
				
				sumX2 = sumX2 + denominator_array(2, (k +199 + helper))
				sumY2 = sumY2 + denominator_array(3, (k  + helper))

			END DO
			Crosscor_factor(i) = (SQRT(sumX * sumY)) / (SQRT(sumX2 * sumY2))

			! Adjust variables
			sumX = 0
			sumX2 = 0
			sumY = 0 
			sumY2 = 0 

	END DO


!	ALLOCATE (Crosscor_factor_helper(segment_numbers, Correlation_steps_mid	))
!	Crosscor_factor_helper = 0

!	DO i=1, segment_numbers
!		helper = segment_rows * (i-1)
!		 DO j=1, Correlation_steps_mid	
!			 DO k = 1, boundary
!				sumX = sumX + (Filtered_array_mid_cut(2 , (k + 199 + helper))) * (Filtered_array_mid_cut(2 , (k + 199 + helper)))
!				sumY = sumY + (Filtered_array_mid_cut(3, (k+ helper+offset))) * (Filtered_array_mid_cut(3, (k+ helper+offset)))
				
!				sumX2 = sumX2 + denominator_array(2, (k +199 + helper))
!				sumY2 = sumY2 + denominator_array(3, (k + helper + offset))

!			END DO 
!			Crosscor_factor_helper(i,j) = (SQRT(sumX * sumY)) / (SQRT(sumX2 * sumY2))

			! Adjust variables
!			offset=offset +1
!			sumX = 0
!			sumX2 = 0
!			sumY = 0 
!			sumY2 = 0 

!		END DO 

!		offset = 0		!zero offset helper
!	END DO 

!	sumX = 0
!	sumX2 = 0
!	sumY = 0 
!	sumY2 = 0 
!	offset = 0
!	helper = 0


	! auto-correlation -> factor
	DO i=1, segment_numbers
		helper = segment_rows * (i-1)
			DO k = 1, boundary
				sumX = sumX + (Filtered_array_mid_cut(2 , (k + helper))) * (Filtered_array_mid_cut(2 , (k + helper)))
	!			sumY = sumY + (Filtered_array_mid_cut(2, (k + helper))) * (Filtered_array_mid_cut(2, (k + helper)))
				
				sumX2 = sumX2 + denominator_array(2, (k + helper))
	!			sumY2 = sumY2 + denominator_array(2, (k  + helper))

			END DO
			Autocor_factor(i) = sumX / sumX2

			! Adjust variables
			sumX = 0
			sumX2 = 0
			sumY = 0 
			sumY2 = 0 
	END DO



!	ALLOCATE (Autocor_factor_helper(segment_numbers, Correlation_steps_mid ))
!	Autocor_factor_helper = 0

!	DO i=1, segment_numbers
!		helper = segment_rows * (i-1)
!		DO j=1, Correlation_steps_mid	
!			DO k = 1, boundary
!				sumX = sumX + (Filtered_array_mid_cut(2 , (k + helper))) * (Filtered_array_mid_cut(2 , (k + helper)))
!				sumY = sumY + (Filtered_array_mid_cut(2, (k+ helper+offset))) * (Filtered_array_mid_cut(2, (k+ helper+offset)))
				
!				sumX2 = sumX2 + denominator_array(2, (k + helper))
!				sumY2 = sumY2 + denominator_array(2, (k + helper + offset))

!			END DO 
!			Autocor_factor_helper(i,j) = (SQRT(sumX * sumY)) / (SQRT(sumX2 * sumY2))

			! Adjust variables
!			offset=offset +1
!			sumX = 0
!			sumX2 = 0
!			sumY = 0 
!			sumY2 = 0 

!		END DO 

!		offset = 0		!zero offset helper
!	END DO 


!	sumX = 0
!	sumX2 = 0
!	sumY = 0 
!	sumY2 = 0 
!	offset = 0
!	helper = 0


	!Performing Auto-correlation
	ALLOCATE (Autocor_helper(segment_numbers, Correlation_steps))
	Autocor_helper = 0			
	outer : DO i=1, segment_numbers
			helper = segment_rows * (i-1)
		middle : DO j=1, Correlation_steps
			inner : DO k = 1, boundary
					sumX = sumX + Filtered_array_mid_cut(2 , (k + helper))
					sumX2 = sumX2 + Filtered_array_mid_cut(2 , (k+ helper))* &
							Filtered_array_mid_cut(2 , (k+ helper)) 
					sumY = sumY + Filtered_array_mid_cut(2, (k+ helper+offset))
					sumY2 = sumY2 + Filtered_array_mid_cut(2, (k+ helper+offset)) &
							*Filtered_array_mid_cut(2, (k+ helper+offset))
					sumXY = sumXY + Filtered_array_mid_cut(2,(k+ helper)) * &
							Filtered_array_mid_cut(2,(k+ helper+offset))
			END DO inner
		
			Autocor_helper(i,j) = (sumXY - sumX * sumY/boundary) /(SQRT(sumX2-sumX*sumX/ &
								  boundary) * SQRT(sumY2-sumY*sumY/boundary))
			! Adjust variables
			offset=offset +1
			sumX = 0
			sumX2 = 0
			sumY = 0 
			sumY2 = 0 
			sumXY = 0 			

		END DO middle
		offset = 0		!zero offset helper
	END DO outer
	helper = 0

	!average the autocorrelation-segments and add time
	DO i = 1, Correlation_steps
			Autocorrelation_mid (1, i) = real((i-1))/frequency
			
			DO j = 1, segment_numbers
					average_value = average_value + (Autocor_helper (j , i)* Autocor_factor(j))
			END DO 
				
				Autocorrelation_mid (2, i) = average_value/segment_numbers
				average_value = 0		
	END DO

	DO j = 1, segment_numbers
	average_factor = average_factor + Autocor_factor(j)
	END Do
	Auto_fac = average_factor/segment_numbers

	average_factor = 0
	

	WRITE (*,*) 'Finished calculation of autocorrelation_mid! '

!		OPEN (UNIT = 2002, FILE = 'auto.txt', STATUS = 'REPLACE', ACTION ='WRITE', IOSTAT=status)
!		WRITE (2002, 2003) ((Autocor_helper (i,j), i=1,segment_numbers), j=1, Correlation_steps)
!		2003 FORMAT (15F10.5)
			
	DEALLOCATE(Autocor_helper)
!	DEALLOCATE(Autocor_factor_helper)
	!Performing Crosscorrelation
	sumX = 0
	sumX2 = 0
	sumY = 0 
	sumY2 = 0 
	sumXY = 0
	offset = 0
	helper = 0


	ALLOCATE (Crosscor_helper(segment_numbers, Correlation_steps))
	Crosscor_helper = 0

	outer_cross : DO i=1, segment_numbers
		helper = segment_rows * (i-1)
		middle_cross : DO j=1, Correlation_steps
			inner_cross : DO k = 1, boundary
				sumX = sumX + Filtered_array_mid_cut(2 , (k + 199 + helper))
				sumX2 = sumX2 + Filtered_array_mid_cut(2 , (k +199 + helper))* &
						Filtered_array_mid_cut(2 , (k +199 + helper)) 
				sumY = sumY + Filtered_array_mid_cut(3, (k+ helper+offset))
				sumY2 = sumY2 + Filtered_array_mid_cut(3, (k+ helper+offset))* &
						Filtered_array_mid_cut(3, (k+ helper+offset))
				sumXY = sumXY + Filtered_array_mid_cut(2, (k + 199 + helper)) * &
						Filtered_array_mid_cut(3, (k+ helper+offset))

			END DO inner_cross
		
			Crosscor_helper(i,j) = (boundary*sumXY - sumX * sumY)/(SQRT(boundary*sumX2- &
									sumX*sumX)* SQRT(boundary*sumY2-sumY*sumY))
			! Adjust variables
			offset=offset +1
			sumX = 0
			sumX2 = 0
			sumY = 0 
			sumY2 = 0 
			sumXY = 0 			

		END DO middle_cross

		offset = 0		!zero offset helper
	END DO outer_cross

	!average the autocorrelation-segments and add time
	DO i = 1, Correlation_steps
			Crosscorrelation_mid (1, i) = -REAL(199)/frequency + real((i-1))/frequency
			
			DO j = 1, segment_numbers
				average_value = average_value + (Crosscor_helper (j , i)* Crosscor_factor(j))
			END DO
					
			Crosscorrelation_mid (2, i) = average_value/segment_numbers
			average_value = 0
	END DO

	DO j = 1, segment_numbers
	average_factor = average_factor + Crosscor_factor(j)
	END Do
	Cross_fac = average_factor/segment_numbers
	 
		WRITE (*,*) 'Finished calculation of crosscorrelation_mid! '

!		OPEN (UNIT = 2000, FILE = 'cross.txt', STATUS = 'REPLACE', ACTION ='WRITE', IOSTAT=status)
!		WRITE (2000, 2001) ((Crosscor_helper (i,j), i=1,segment_numbers), j=1, Correlation_steps)
!		2001 FORMAT (15F10.5)
			
	DEALLOCATE(Crosscor_helper)
!	DEALLOCATE(Crosscor_factor_helper)
	DEALLOCATE(denominator_array)
	
		
END SUBROUTINE Correlation_mid


!_______________________________________________________________________________________________________________


SUBROUTINE Correlation_scales_mid (Autocorrelation_mid, Crosscorrelation_mid, Correlation_steps, &
		   CorScale_results_mid, frequency, delta_x, velocity, Auto_fac, Cross_fac)
IMPLICIT NONE

	INTEGER, PARAMETER :: DBL = SELECTED_REAL_KIND(p=15)	!Precision of the Real variables
	INTEGER, INTENT (IN) :: Correlation_steps
	INTEGER, INTENT (IN) :: frequency
	REAL, INTENT (IN) :: delta_x
	REAL(kind=DBL), INTENT (INOUT) :: velocity
	REAL(kind=DBL), INTENT (INOUT) :: Auto_fac
	REAL(kind=DBL), INTENT (INOUT) :: Cross_fac
	REAL(kind=DBL), INTENT (IN) :: Autocorrelation_mid (2, Correlation_steps)
	REAL(kind=DBL), INTENT (IN) :: Crosscorrelation_mid (2, Correlation_steps)
	REAL(kind=DBL), INTENT (OUT) ::CorScale_results_mid(18)
	INTEGER :: i = 0 , j = 0, counter, locator
	REAL(kind=DBL) :: helper
	INTEGER :: helper_int
	REAL(kind=DBL) :: min_auto, min_auto_loc, min_cross, min_cross_loc 
	REAL(kind=DBL) :: T_point_five
	REAL(kind=DBL) :: Zero_crossing_auto
	REAL(kind=DBL) :: Txx
	REAL(kind=DBL) :: Rxz_max
	REAL(kind=DBL) :: Rxx_max
	REAL(kind=DBL) :: helper_RxxMax_half
	REAL(kind=DBL) :: T_Rxz_max
	REAL(kind=DBL) :: T_half_Rxz_max 
	REAL(kind=DBL) :: Zero_crossing_cross 
	REAL(kind=DBL) :: Txz
	REAL(kind=DBL) :: Tu
	REAL(kind=DBL) :: v_prime 

	velocity = 0
	CorScale_results_mid = 0
	helper = 0
	T_point_five = 0
	Zero_crossing_auto =0
	min_auto_loc = 0
	min_auto = 0
	Txx = 0
	Rxz_max = 0
	T_Rxz_max = 0
	Zero_crossing_cross = 0
	min_cross_loc = 0
	min_cross = 0
	T_half_Rxz_max = 0
	Txz = 0
	Tu = 0
	Velocity = 0
	counter=0
	locator = 0
	helper_int = 0
	v_prime = 0
	Rxx_max = 0

	
	Rxx_max  = Autocorrelation_mid (2,1)
	helper_RxxMax_half = Rxx_max/2

	!Autocorrelation_results
	!Calculation of Time when Autocorrelation is 0.5
	DO i=1, Correlation_steps
		helper = Autocorrelation_mid(2, i)
		IF (helper == helper_RxxMax_half) THEN
			T_point_five = Autocorrelation_mid(1,i)
		EXIT
			ELSE IF (helper < helper_RxxMax_half) THEN
			    T_point_five = Autocorrelation_mid(1,i-1) + (Autocorrelation_mid(1,i) - &
				Autocorrelation_mid(1,i-1)) / (Autocorrelation_mid(2,i) - Autocorrelation_mid(2,i-1))&
				 * (helper_RxxMax_half-Autocorrelation_mid(2,i-1))	
			EXIT
		END IF
	END DO

	helper=0
	counter = 0

	!Calculation of Zero crossing of Autocorrelation function
	DO i=1, Correlation_steps
		helper = Autocorrelation_mid(2,i)
		Zero_crossing_auto = 0
		counter = counter +1
		IF (helper == 0) THEN
			Zero_crossing_auto = Autocorrelation_mid(1,i)
				EXIT
			ELSE IF (helper < 0.) THEN
				Zero_crossing_auto = Autocorrelation_mid(1,i-1) + (Autocorrelation_mid(1,i) - &
				Autocorrelation_mid(1,i-1)) / (Autocorrelation_mid(2,i) - Autocorrelation_mid(2,i-1)) &
				* (0-Autocorrelation_mid(2,i-1))
				EXIT
				ELSE 
				Zero_crossing_auto = 99	
		END IF
	END DO	

	!Calculation of Txx
	min_auto=1
	locator = 0
	min_auto_loc = 0
	Txx =0
	IF (INT(Zero_crossing_auto) == 99) THEN
		DO j = 1, Correlation_steps-1
			min_auto = Min(Autocorrelation_mid(2,j), min_auto)
		END DO

		DO i = 1, Correlation_steps-1
		locator = locator +1
			IF (Autocorrelation_mid(2,i) - min_auto < 0.00001) THEN
			min_auto_loc = Autocorrelation_mid(1,locator)
			EXIT
			END IF
		END DO

		DO i = 1, locator
			IF (i ==1000) THEN
			Txx = Txx + (Autocorrelation_mid(2,i)/frequency)
			ELSE
			Txx = Txx + (Autocorrelation_mid(2,i)+Autocorrelation_mid(2,i+1))/(2*frequency)
			END IF
		END DO
	ELSE 
		DO j = 1, counter
			IF (j ==1000) THEN
			Txx = Txx + (Autocorrelation_mid(2,j)/frequency)
			ELSE
			Txx = Txx + (Autocorrelation_mid(2,j)+Autocorrelation_mid(2,j+1))/(2*frequency)
			END IF
		END DO 
		min_auto=0
		min_auto_loc=Zero_crossing_auto
	END IF
		
	helper =0
	counter =0
	Rxz_max = 0
	!Crosscorrelation_results
	!Find maximum Crosscorrelation value (Rxz)max
	DO i = 1, Correlation_steps
		helper = Crosscorrelation_mid(2,i)
		Rxz_max = Max(Rxz_max,helper)
	END DO	
		
	helper = 0 
	locator = 0
	T_Rxz_max = 0

	!Corresponding Time for Rxz_max
	DO i = 1, Correlation_steps
		locator = locator +1
		IF (Rxz_max - Crosscorrelation_mid(2,i) < 0.0000001) THEN
		T_Rxz_max = Crosscorrelation_mid(1,locator)
		EXIT 
		END IF
	END DO

	T_half_Rxz_max = 0

	!Calculation of Time for 0.5*Rxz_max
	DO i=locator, Correlation_steps
		IF (Crosscorrelation_mid (2,i) <= 0.5*Rxz_max) THEN
			T_half_Rxz_max = Crosscorrelation_mid(1,i-1) + (Crosscorrelation_mid(1,i) - &
							Crosscorrelation_mid(1,i-1)) / (Crosscorrelation_mid(2,i) - &
							Crosscorrelation_mid(2,i-1)) * (0.5*Rxz_max-Crosscorrelation_mid(2,i-1))
			EXIT
		END IF
	END DO

	counter =-1
	helper = 0

	!Time where the values of Rxz_max cross the x-axis for the first time right of Rxy_max
	DO i=locator, Correlation_steps
		helper = Crosscorrelation_mid(2,i)
		Zero_crossing_cross = 0
		counter = counter +1
		IF (helper == 0) THEN
			Zero_crossing_cross = Crosscorrelation_mid(1,i)
				EXIT
			ELSE IF (helper < 0.) THEN
				Zero_crossing_cross = Crosscorrelation_mid(1,i-1) + (Crosscorrelation_mid(1,i) - &
			    Crosscorrelation_mid(1,i-1))  / (Crosscorrelation_mid(2,i) - Crosscorrelation_mid(2,i-1))&
				* (0-Crosscorrelation_mid(2,i-1))
				EXIT
				ELSE 
				Zero_crossing_cross = 98	
		END IF
	END DO	

	!Calculation of Txz
	min_cross = Rxz_max
	min_cross_loc = 0
	Txz = 0
	IF (INT(Zero_crossing_cross) == 98) THEN
		DO j = locator, Correlation_steps-1
			min_cross = Min(Crosscorrelation_mid(2,j), min_cross)
		END DO
		
		helper_int = 0
		DO i = locator, Correlation_steps-1
		helper_int = helper_int +1
			IF (Crosscorrelation_mid(2,i) - min_cross < 0.0000001) THEN
			min_cross_loc = Crosscorrelation_mid(1,(helper_int+locator-1))
			EXIT
			END IF
		END DO

		DO i = locator, (helper_int +locator-1)
			IF (i ==1000) THEN
			Txz = Txz + (Crosscorrelation_mid(2,i)/frequency)
			ELSE
			Txz = Txz + (Crosscorrelation_mid(2,i)+Crosscorrelation_mid(2,i+1))/(2*frequency)
			END IF
		END DO
	ELSE 
		DO j = locator, (locator+counter)
			IF (j ==1000) THEN
			Txz = Txz + (Crosscorrelation_mid(2,j)/frequency)
			ELSE
			Txz = Txz + (Crosscorrelation_mid(2,j)+Crosscorrelation_mid(2,j+1))/(2*frequency)
			END IF
		END DO 
		min_cross = 0
		min_cross_loc = Zero_crossing_cross
	END IF
		
	!Calculation of turbulence intensity
	IF (((T_half_Rxz_max - T_Rxz_max)*(T_half_Rxz_max - T_Rxz_max)- &
		 T_point_five*T_point_five) <0) THEN
		 Tu = 97
	ELSE
	Tu = 0.851 * SQRT((T_half_Rxz_max - T_Rxz_max)*(T_half_Rxz_max - T_Rxz_max)- &
		 T_point_five*T_point_five)/T_Rxz_max
	END IF
		
	!Calculation of time averaged interfacial velocity
	velocity = delta_x/1000/T_Rxz_max	

	v_prime = Tu * velocity

	!store results into array CorScale_results  :T_point_five, Zero_crossing_auto, min_auto_loc, min_auto,Txx, &
						!T_Rxz_max, T_half_Rxz_max, Zero_crossing_cross, min_cross_loc, min_cross,Txz, Tu, Velocity, Rxx_max
	CorScale_results_mid = (/ T_point_five, Zero_crossing_auto, min_auto_loc, min_auto, Txx, Rxz_max, T_Rxz_max, &
						Zero_crossing_cross, min_cross_loc, min_cross, T_half_Rxz_max, Txz, Tu, Velocity, v_prime, Rxx_max, Auto_fac, Cross_fac /)


END SUBROUTINE Correlation_scales_mid

!____________________________________
!______________________________________________________________________________________________________


SUBROUTINE Correlation_low (Filtered_array_mid_cut, Filtered_array_high_cut, number_columns, number_rows_cut, Correlation_steps_low, & 
			         Autocorrelation_fast_low, Crosscorrelation_fast_low, Autocorrelation_low_fast, Crosscorrelation_low_fast, &
					 frequency, segment_numbers, Factors_crossproduct)
IMPLICIT NONE

	INTEGER, PARAMETER :: DBL = SELECTED_REAL_KIND(p=15)	!Precision of the Real variables
	INTEGER, INTENT (IN) :: number_columns
	INTEGER, INTENT (IN) :: number_rows_cut
	INTEGER, INTENT (IN) :: Correlation_steps_low
	INTEGER, INTENT (IN) :: frequency
	INTEGER, INTENT (IN) :: segment_numbers
	REAL(kind=DBL), INTENT (IN) :: Filtered_array_mid_cut (number_columns, number_rows_cut)
	REAL(kind=DBL), INTENT (IN) :: Filtered_array_high_cut (number_columns, number_rows_cut)
	REAL(kind=DBL), INTENT (INOUT) :: Factors_crossproduct(3)
	REAL(kind=DBL), INTENT (OUT) :: Autocorrelation_fast_low (2, Correlation_steps_low)
	REAL(kind=DBL), INTENT (OUT) :: Crosscorrelation_fast_low(2, Correlation_steps_low)
	REAL(kind=DBL), INTENT (OUT) :: Autocorrelation_low_fast (2, Correlation_steps_low)
	REAL(kind=DBL), INTENT (OUT) :: Crosscorrelation_low_fast(2, Correlation_steps_low)
	INTEGER :: segment_rows							!number of rows in each segment
	INTEGER :: helper								!helper for faster correlations
	INTEGER :: boundary								!helper for faster correlation
	INTEGER :: i =0, j = 0, k=0							!integer for loop
	REAL(kind=DBL), ALLOCATABLE, DIMENSION  (:,:) :: Autocor_helper	
	REAL(kind=DBL), ALLOCATABLE, DIMENSION  (:,:) :: Crosscor_helper	
	REAL(kind=DBL) :: sumX, sumY, sumX2, sumY2, sumXY !Variables for correl
	INTEGER :: offset												 !integer for offsetting
	INTEGER :: status
	REAL(kind=DBL) :: average_value				!Helper to calculate the average values 
	
	
	! calculation for factor
	REAL(kind=DBL), ALLOCATABLE, DIMENSION  (:,:) :: denominator_array	
	REAL(kind=DBL), ALLOCATABLE, DIMENSION  (:,:) :: Crosscor_factor_helper	
	REAL(kind=DBL), ALLOCATABLE, DIMENSION  (:,:) :: Autocor_factor_helper	

	REAL(kind=DBL) :: Crosscor_factor(segment_numbers)
	REAL(kind=DBL) :: Autocor_factor(segment_numbers)
	REAL(kind=DBL) :: Boundary_factor
	REAL(kind=DBL) :: 	average_factor

	segment_rows	=0
	helper = 0
	boundary = 0
	sumX = 0
	sumY = 0
	sumX2 = 0
	sumY2 = 0
	sumXY = 0
	offset = 0
	status = 0
	average_value =0
	Crosscor_factor = 0 
	Autocor_factor = 0 
	Boundary_factor = 0
	average_factor = 0

	segment_rows = number_rows_cut/segment_numbers
	boundary = segment_rows - Correlation_steps_low

	Boundary_factor = number_rows_cut - Correlation_steps_low	


	!FACTOR calculations
	ALLOCATE (denominator_array(number_columns, number_rows_cut))
	denominator_array = 0

		Do j = 1, number_rows_cut
			denominator_array (2,j) =  (Filtered_array_mid_cut (2, j) + Filtered_array_high_cut (2, j)) * (Filtered_array_mid_cut (2, j) + Filtered_array_high_cut (2, j))
			denominator_array (3,j) =  (Filtered_array_mid_cut (3, j) + Filtered_array_high_cut (3, j)) * (Filtered_array_mid_cut (3, j) + Filtered_array_high_cut (3, j))
		END DO


	! auto-correlation -> factor fast_low
	DO i=1, segment_numbers
		helper = segment_rows * (i-1)
			DO k = 1, boundary
				sumX = sumX + (Filtered_array_high_cut(2 , (k + helper))) * (Filtered_array_high_cut(2 , (k + helper)))
				sumY = sumY + (Filtered_array_mid_cut(2, (k + helper))) * (Filtered_array_mid_cut(2, (k + helper)))
				
				sumX2 = sumX2 + denominator_array(2, (k + helper))
	!			sumY2 = sumY2 + denominator_array(2, (k  + helper))

			END DO
			Autocor_factor(i) = (SQRT(sumX * sumY)) / sumX2

			! Adjust variables
			sumX = 0
			sumX2 = 0
			sumY = 0 
			sumY2 = 0 
	END DO


	!Performing Auto-correlation
	ALLOCATE (Autocor_helper(segment_numbers, Correlation_steps_low))
	Autocor_helper = 0			
	outer : DO i=1, segment_numbers
			helper = segment_rows * (i-1)
		middle : DO j=1, Correlation_steps_low
			inner : DO k = 1, boundary
					sumX = sumX + Filtered_array_high_cut(2 , (k + helper))
					sumX2 = sumX2 + Filtered_array_high_cut(2 , (k+ helper))* &
							Filtered_array_high_cut(2 , (k+ helper)) 
					sumY = sumY + Filtered_array_mid_cut(2, (k+ helper+offset))
					sumY2 = sumY2 + Filtered_array_mid_cut(2, (k+ helper+offset)) &
							*Filtered_array_mid_cut(2, (k+ helper+offset))
					sumXY = sumXY + Filtered_array_high_cut(2,(k+ helper)) * &
							Filtered_array_mid_cut(2,(k+ helper+offset))
			END DO inner
		
			Autocor_helper(i,j) = (sumXY - sumX * sumY/boundary) /(SQRT(sumX2-sumX*sumX/ &
								  boundary) * SQRT(sumY2-sumY*sumY/boundary))
			! Adjust variables
			offset=offset +1
			sumX = 0
			sumX2 = 0
			sumY = 0 
			sumY2 = 0 
			sumXY = 0 			

		END DO middle
		offset = 0		!zero offset helper
	END DO outer
	helper = 0

	!average the autocorrelation-segments and add time
	DO i = 1, Correlation_steps_low
			Autocorrelation_fast_low (1, i) = real((i-1))/frequency
			
			DO j = 1, segment_numbers
					average_value = average_value + (Autocor_helper (j , i)* Autocor_factor(j))
			END DO 
				
				Autocorrelation_fast_low (2, i) = average_value/segment_numbers
				average_value = 0
				
	END DO


	WRITE (*,*) 'Finished calculation of Autocorrelation_fast_low! '

!		OPEN (UNIT = 2002, FILE = 'auto.txt', STATUS = 'REPLACE', ACTION ='WRITE', IOSTAT=status)
!		WRITE (2002, 2003) ((Autocor_helper (i,j), i=1,segment_numbers), j=1, Correlation_steps)
!		2003 FORMAT (15F10.5)


	sumX = 0
	sumX2 = 0
	sumY = 0 
	sumY2 = 0 
	sumXY = 0
	offset = 0
	helper = 0
	Autocor_helper = 0	
	average_value = 0

	!Performing Auto-correlation for slow_fast
	DO i=1, segment_numbers
			helper = segment_rows * (i-1)
		DO j=1, Correlation_steps_low
			DO k = 1, boundary
					sumX = sumX + Filtered_array_mid_cut(2 , (k + helper))
					sumX2 = sumX2 + Filtered_array_mid_cut(2 , (k+ helper))* &
							Filtered_array_mid_cut(2 , (k+ helper)) 
					sumY = sumY + Filtered_array_high_cut(2, (k+ helper+offset))
					sumY2 = sumY2 + Filtered_array_high_cut(2, (k+ helper+offset)) &
							*Filtered_array_high_cut(2, (k+ helper+offset))
					sumXY = sumXY + Filtered_array_mid_cut(2,(k+ helper)) * &
							Filtered_array_high_cut(2,(k+ helper+offset))
			END DO
		
			Autocor_helper(i,j) = (sumXY - sumX * sumY/boundary) /(SQRT(sumX2-sumX*sumX/ &
								  boundary) * SQRT(sumY2-sumY*sumY/boundary))
			! Adjust variables
			offset=offset +1
			sumX = 0
			sumX2 = 0
			sumY = 0 
			sumY2 = 0 
			sumXY = 0 			

		END DO
		offset = 0		!zero offset helper
	END DO
	helper = 0

	!average the autocorrelation-segments and add time
	DO i = 1, Correlation_steps_low
			Autocorrelation_low_fast (1, i) = real((i-1))/frequency
			
			DO j = 1, segment_numbers
					average_value = average_value + (Autocor_helper (j , i)* Autocor_factor(j))
			END DO 
				
				Autocorrelation_low_fast (2, i) = average_value/segment_numbers
				average_value = 0
				
	END DO

	DO j = 1, segment_numbers
	average_factor = average_factor + Autocor_factor(j)
	END DO
	Factors_crossproduct(1) = average_factor/segment_numbers

	average_factor = 0

	WRITE (*,*) 'Finished calculation of Autocorrelation_low_fast! '


	DEALLOCATE(Autocor_helper)
!	DEALLOCATE(Autocor_factor_helper)
	!Performing Crosscorrelation
	sumX = 0
	sumX2 = 0
	sumY = 0 
	sumY2 = 0 
	sumXY = 0
	offset = 0
	helper = 0
	average_value = 0



	! cross-correlation -> factor fast_low
	DO i=1, segment_numbers
		helper = segment_rows * (i-1)
			DO k = 1, boundary
				sumX = sumX + (Filtered_array_high_cut(2 , (k + 199 + helper))) * (Filtered_array_high_cut(2 , (k + 199 + helper)))
				sumY = sumY + (Filtered_array_mid_cut(3, (k + helper))) * (Filtered_array_mid_cut(3, (k + helper)))
				
				sumX2 = sumX2 + denominator_array(2, (k +199 + helper))
				sumY2 = sumY2 + denominator_array(3, (k  + helper))
			END DO
			Crosscor_factor(i) = (SQRT(sumX * sumY)) / (SQRT(sumX2 * sumY2))

			! Adjust variables
			sumX = 0
			sumX2 = 0
			sumY = 0 
			sumY2 = 0 

	END DO

	ALLOCATE (Crosscor_helper(segment_numbers, Correlation_steps_low))
	Crosscor_helper = 0

	DO i=1, segment_numbers
		helper = segment_rows * (i-1)
		DO j=1, Correlation_steps_low
			DO k = 1, boundary
				sumX = sumX + Filtered_array_high_cut(2 , (k + 199 + helper))
				sumX2 = sumX2 + Filtered_array_high_cut(2 , (k +199 + helper))* &
						Filtered_array_high_cut(2 , (k +199 + helper)) 
				sumY = sumY + Filtered_array_mid_cut(3, (k+ helper+offset))
				sumY2 = sumY2 + Filtered_array_mid_cut(3, (k+ helper+offset))* &
						Filtered_array_mid_cut(3, (k+ helper+offset))
				sumXY = sumXY + Filtered_array_high_cut(2, (k + 199 + helper)) * &
						Filtered_array_mid_cut(3, (k+ helper+offset))

			END DO 
		
			Crosscor_helper(i,j) = (boundary*sumXY - sumX * sumY)/(SQRT(boundary*sumX2- &
									sumX*sumX)* SQRT(boundary*sumY2-sumY*sumY))
			! Adjust variables
			offset=offset +1
			sumX = 0
			sumX2 = 0
			sumY = 0 
			sumY2 = 0 
			sumXY = 0 			

		END DO 
		offset = 0		!zero offset helper
	END DO

	!average the crosscorrelation-segments and add time
	DO i = 1, Correlation_steps_low
			Crosscorrelation_fast_low (1, i) = -REAL(199)/frequency + real((i-1))/frequency
			
			DO j = 1, segment_numbers
				average_value = average_value + (Crosscor_helper (j , i)* Crosscor_factor(j))
			END DO
					
			Crosscorrelation_fast_low (2, i) = average_value/segment_numbers
			average_value = 0
	END DO

	DO j = 1, segment_numbers
	average_factor = average_factor + Crosscor_factor(j)
	END DO
	Factors_crossproduct(3) = average_factor/segment_numbers

	average_factor = 0

	 
		WRITE (*,*) 'Finished calculation of Crosscorrelation_fast_low! '

!		OPEN (UNIT = 2000, FILE = 'cross.txt', STATUS = 'REPLACE', ACTION ='WRITE', IOSTAT=status)
!		WRITE (2000, 2001) ((Crosscor_helper (i,j), i=1,segment_numbers), j=1, Correlation_steps)
!		2001 FORMAT (15F10.5)

	Crosscor_factor = 0
	Crosscor_helper = 0
	average_value = 0

	! cross-correlation -> factor low_fast
	DO i=1, segment_numbers
		helper = segment_rows * (i-1)
			DO k = 1, boundary
				sumX = sumX + (Filtered_array_mid_cut(2 , (k + 199 + helper))) * (Filtered_array_mid_cut(2 , (k + 199 + helper)))
				sumY = sumY + (Filtered_array_high_cut(3, (k + helper))) * (Filtered_array_high_cut(3, (k + helper)))
				
				sumX2 = sumX2 + denominator_array(2, (k +199 + helper))
				sumY2 = sumY2 + denominator_array(3, (k  + helper))
			END DO
			Crosscor_factor(i) = (SQRT(sumX * sumY)) / (SQRT(sumX2 * sumY2))

			! Adjust variables
			sumX = 0
			sumX2 = 0
			sumY = 0 
			sumY2 = 0 

	END DO

	DO i=1, segment_numbers
		helper = segment_rows * (i-1)
		DO j=1, Correlation_steps_low
			DO k = 1, boundary
				sumX = sumX + Filtered_array_mid_cut(2 , (k + 199 + helper))
				sumX2 = sumX2 + Filtered_array_mid_cut(2 , (k +199 + helper))* &
						Filtered_array_mid_cut(2 , (k +199 + helper)) 
				sumY = sumY + Filtered_array_high_cut(3, (k+ helper+offset))
				sumY2 = sumY2 + Filtered_array_high_cut(3, (k+ helper+offset))* &
						Filtered_array_high_cut(3, (k+ helper+offset))
				sumXY = sumXY + Filtered_array_mid_cut(2, (k + 199 + helper)) * &
						Filtered_array_high_cut(3, (k+ helper+offset))

			END DO 
		
			Crosscor_helper(i,j) = (boundary*sumXY - sumX * sumY)/(SQRT(boundary*sumX2- &
									sumX*sumX)* SQRT(boundary*sumY2-sumY*sumY))
			! Adjust variables
			offset=offset +1
			sumX = 0
			sumX2 = 0
			sumY = 0 
			sumY2 = 0 
			sumXY = 0 			

		END DO 
		offset = 0		!zero offset helper
	END DO

	!average the crosscorrelation-segments and add time
	DO i = 1, Correlation_steps_low
			Crosscorrelation_low_fast (1, i) = -REAL(199)/frequency + real((i-1))/frequency
			
			DO j = 1, segment_numbers
				average_value = average_value + (Crosscor_helper (j , i)* Crosscor_factor(j))
			END DO
					
			Crosscorrelation_low_fast (2, i) = average_value/segment_numbers
			average_value = 0
	END DO

	DO j = 1, segment_numbers
	average_factor = average_factor + Crosscor_factor(j)
	END DO
	Factors_crossproduct(2) = average_factor/segment_numbers

	WRITE (*,*) 'Finished calculation of Crosscorrelation_low_fast! '
	
	DEALLOCATE(Crosscor_helper)
!	DEALLOCATE(Crosscor_factor_helper)
	DEALLOCATE(denominator_array)

		
END SUBROUTINE Correlation_low

!_______________________________________________________________________________________________________________


SUBROUTINE Correlation_scales_low (Autocorrelation_low, Crosscorrelation_low, Correlation_steps, &
		   CorScale_results_low, frequency, delta_x, velocity)
IMPLICIT NONE

	INTEGER, PARAMETER :: DBL = SELECTED_REAL_KIND(p=15)	!Precision of the Real variables
	INTEGER, INTENT (IN) :: Correlation_steps
	INTEGER, INTENT (IN) :: frequency
	REAL, INTENT (IN) :: delta_x
	REAL(kind=DBL), INTENT (INOUT) :: velocity
	REAL(kind=DBL), INTENT (IN) :: Autocorrelation_low (2, Correlation_steps)
	REAL(kind=DBL), INTENT (IN) :: Crosscorrelation_low (2, Correlation_steps)
	REAL(kind=DBL), INTENT (OUT) ::CorScale_results_low(16)
	INTEGER :: i = 0 , j = 0, counter, locator
	REAL(kind=DBL) :: helper
	INTEGER :: helper_int
	REAL(kind=DBL) :: min_auto, min_auto_loc, min_cross, min_cross_loc 
	REAL(kind=DBL) :: T_point_five
	REAL(kind=DBL) :: Zero_crossing_auto
	REAL(kind=DBL) :: Txx
	REAL(kind=DBL) :: Rxz_max
	REAL(kind=DBL) :: Rxx_max
	REAL(kind=DBL) :: T_Rxz_max
	REAL(kind=DBL) :: T_half_Rxz_max 
	REAL(kind=DBL) :: Zero_crossing_cross 
	REAL(kind=DBL) :: Txz
	REAL(kind=DBL) :: Tu 
	REAL(kind=DBL) :: v_prime
	REAL(kind=DBL) :: helper_RxxMax_half

	velocity = 0
	CorScale_results_low = 0
	helper = 0
	T_point_five = 0
	Zero_crossing_auto =0
	min_auto_loc = 0
	min_auto = 0
	Txx = 0
	Rxz_max = 0
	T_Rxz_max = 0
	Zero_crossing_cross = 0
	min_cross_loc = 0
	min_cross = 0
	T_half_Rxz_max = 0
	Txz = 0
	Tu = 0
	Velocity = 0
	counter=0
	locator = 0
	helper_int = 0
	v_prime = 0
	Rxx_max = 0
	
	Rxx_max  = Autocorrelation_low (2,1)
	helper_RxxMax_half = Rxx_max/2

	!Autocorrelation_results
	!Calculation of Time when Autocorrelation is 0.5
	DO i=1, Correlation_steps
		helper = Autocorrelation_low(2, i)
		IF (helper == helper_RxxMax_half) THEN
			T_point_five = Autocorrelation_low(1,i)
		EXIT
			ELSE IF (helper < helper_RxxMax_half) THEN
			    T_point_five = Autocorrelation_low(1,i-1) + (Autocorrelation_low(1,i) - &
				Autocorrelation_low(1,i-1)) / (Autocorrelation_low(2,i) - Autocorrelation_low(2,i-1))&
				 * (helper_RxxMax_half - Autocorrelation_low(2,i-1))	
			EXIT
		END IF
	END DO

	helper=0
	counter = 0

	!Calculation of Zero crossing of Autocorrelation function
	DO i=1, Correlation_steps
		helper = Autocorrelation_low(2,i)
		Zero_crossing_auto = 0
		counter = counter +1
		IF (helper == 0) THEN
			Zero_crossing_auto = Autocorrelation_low(1,i)
				EXIT
			ELSE IF (helper < 0.) THEN
				Zero_crossing_auto = Autocorrelation_low(1,i-1) + (Autocorrelation_low(1,i) - &
				Autocorrelation_low(1,i-1)) / (Autocorrelation_low(2,i) - Autocorrelation_low(2,i-1)) &
				* (0-Autocorrelation_low(2,i-1))
				EXIT
				ELSE 
				Zero_crossing_auto = 99	
		END IF
	END DO	

	!Calculation of Txx
	min_auto=1
	locator = 0
	min_auto_loc = 0
	Txx =0
	IF (INT(Zero_crossing_auto) == 99) THEN
		DO j = 1, Correlation_steps-1
			min_auto = Min(Autocorrelation_low(2,j), min_auto)
		END DO

		DO i = 1, Correlation_steps-1
		locator = locator +1
			IF (Autocorrelation_low(2,i) - min_auto < 0.0000001) THEN
			min_auto_loc = Autocorrelation_low(1,locator)
			EXIT
			END IF
		END DO

		DO i = 1, locator
			IF (i ==1000) THEN
			Txx = Txx + (Autocorrelation_low(2,i)/frequency)
			ELSE
			Txx = Txx + (Autocorrelation_low(2,i)+Autocorrelation_low(2,i+1))/(2*frequency)
			END IF
		END DO
	ELSE 
		DO j = 1, counter
			IF (j ==1000) THEN
			Txx = Txx + (Autocorrelation_low(2,j)/frequency)
			ELSE
			Txx = Txx + (Autocorrelation_low(2,j)+Autocorrelation_low(2,j+1))/(2*frequency)
			END IF
		END DO 
		min_auto=0
		min_auto_loc=Zero_crossing_auto
	END IF
		
	helper =0
	counter =0
	Rxz_max = 0
	!Crosscorrelation_results
	!Find maximum Crosscorrelation value (Rxz)max
	DO i = 1, Correlation_steps
		helper = Crosscorrelation_low(2,i)
		Rxz_max = Max(Rxz_max,helper)
	END DO	
		
	helper = 0 
	locator = 0
	T_Rxz_max = 0

	!Corresponding Time for Rxz_max
	DO i = 1, Correlation_steps
		locator = locator +1
		IF (Rxz_max - Crosscorrelation_low(2,i) < 0.0000000001) THEN
		T_Rxz_max = Crosscorrelation_low(1,locator)
		EXIT 
		END IF
	END DO

	T_half_Rxz_max = 0

	!Calculation of Time for 0.5*Rxz_max
	DO i=locator, Correlation_steps
		IF (Crosscorrelation_low (2,i) <= 0.5*Rxz_max) THEN
			T_half_Rxz_max = Crosscorrelation_low(1,i-1) + (Crosscorrelation_low(1,i) - &
							Crosscorrelation_low(1,i-1)) / (Crosscorrelation_low(2,i) - &
							Crosscorrelation_low(2,i-1)) * (0.5*Rxz_max-Crosscorrelation_low(2,i-1))
			EXIT
		END IF
	END DO

	counter =-1
	helper = 0

	!Time where the values of Rxz_max cross the x-axis for the first time right of Rxy_max
	DO i=locator, Correlation_steps
		helper = Crosscorrelation_low(2,i)
		Zero_crossing_cross = 0
		counter = counter +1
		IF (helper == 0) THEN
			Zero_crossing_cross = Crosscorrelation_low(1,i)
				EXIT
			ELSE IF (helper < 0.) THEN
				Zero_crossing_cross = Crosscorrelation_low(1,i-1) + (Crosscorrelation_low(1,i) - &
			    Crosscorrelation_low(1,i-1))  / (Crosscorrelation_low(2,i) - Crosscorrelation_low(2,i-1))&
				* (0-Crosscorrelation_low(2,i-1))
				EXIT
				ELSE 
				Zero_crossing_cross = 98	
		END IF
	END DO	

	!Calculation of Txz
	min_cross = Rxz_max
	min_cross_loc = 0
	Txz = 0
	IF (INT(Zero_crossing_cross) == 98) THEN
		DO j = locator, Correlation_steps-1
			min_cross = Min(Crosscorrelation_low(2,j), min_cross)
		END DO
		
		helper_int = 0
		DO i = locator, Correlation_steps-1
		helper_int = helper_int +1
			IF (Crosscorrelation_low(2,i) - min_cross < 0.0000001) THEN
			min_cross_loc = Crosscorrelation_low(1,(helper_int+locator-1))
			EXIT
			END IF
		END DO

		DO i = locator, (helper_int +locator-1)
			IF (i ==1000) THEN
			Txz = Txz + (Crosscorrelation_low(2,i)/frequency)
			ELSE
			Txz = Txz + (Crosscorrelation_low(2,i)+Crosscorrelation_low(2,i+1))/(2*frequency)
			END IF
		END DO
	ELSE 
		DO j = locator, (locator+counter)
			IF (j ==1000) THEN
			Txz = Txz + (Crosscorrelation_low(2,j)/frequency)
			ELSE
			Txz = Txz + (Crosscorrelation_low(2,j)+Crosscorrelation_low(2,j+1))/(2*frequency)
			END IF
		END DO 
		min_cross = 0
		min_cross_loc = Zero_crossing_cross
	END IF
		
	!Calculation of turbulence intensity
	IF (((T_half_Rxz_max - T_Rxz_max)*(T_half_Rxz_max - T_Rxz_max)- &
		 T_point_five*T_point_five) <0) THEN
		 Tu = 97
	ELSE
	Tu = 0.851 * SQRT((T_half_Rxz_max - T_Rxz_max)*(T_half_Rxz_max - T_Rxz_max)- &
		 T_point_five*T_point_five)/T_Rxz_max
	END IF
		
	!Calculation of time averaged interfacial velocity
	velocity = delta_x/1000/T_Rxz_max
	
	v_prime = Tu * velocity	

	!store results into array CorScale_results  :T_point_five, Zero_crossing_auto, min_auto_loc, min_auto,Txx, &
						!T_Rxz_max, T_half_Rxz_max, Zero_crossing_cross, min_cross_loc, min_cross,Txz, Tu, Velocity,Rxx_max
	CorScale_results_low = (/ T_point_five, Zero_crossing_auto, min_auto_loc, min_auto, Txx, Rxz_max, T_Rxz_max, &
						Zero_crossing_cross, min_cross_loc, min_cross, T_half_Rxz_max, Txz, Tu, Velocity, v_prime, Rxx_max /)

END SUBROUTINE Correlation_scales_low

!____________________________________

SUBROUTINE Probability_V (Data_array, number_columns, number_rows, PDF_V_segment_size, & 
						  PDF_segments, PDF_V, threshold, threshold_C) 

IMPLICIT NONE

INTEGER, PARAMETER :: DBL = SELECTED_REAL_KIND(p=15)	!Precision of the Real variables
	INTEGER, INTENT (IN) :: number_columns
	INTEGER, INTENT (IN) :: number_rows
	REAL, INTENT (IN) :: PDF_V_segment_size
	REAL(kind=DBL), INTENT (IN) :: PDF_segments
	REAL, INTENT (IN) :: threshold_C
	REAL(kind=DBL), INTENT (IN) :: Data_array (number_columns, number_rows)
	REAL(kind=DBL), INTENT (INOUT) :: PDF_V (number_columns, Int(PDF_segments))
	REAL(kind=DBL), INTENT (OUT) :: threshold(number_columns-1)
	INTEGER, ALLOCATABLE, DIMENSION  (:,:) :: PDF_helper
	INTEGER :: i, j, k, l, locator
	REAL(kind=DBL) :: segmenter
	REAL(kind=DBL) :: Max_air, Max_water
	REAL(kind=DBL) :: max_air_loc, max_water_loc
	INTEGER :: PDF_upper_boundary_air

	threshold = 0
	locator = 0
	segmenter = 0
	Max_air = 0
	Max_water = 0
	max_air_loc = 0
	max_water_loc = 0
	PDF_upper_boundary_air = 0

	PDF_upper_boundary_air = 4.0 / PDF_V_segment_size

	ALLOCATE (PDF_helper(number_columns-1, Int(PDF_segments)))
	PDF_helper = 0

	DO i = 1, (number_columns-1)				!number of devices
		DO j = 1, number_rows					!number of rows in data file
			DO k = 1, Int(PDF_segments)			!
				segmenter = k * PDF_V_segment_size
				IF (Data_array(i+1, j) < segmenter) THEN 
					PDF_helper(i,k) = PDF_helper(i,k) +1
					EXIT
				END IF
			END DO
		END DO
	END DO

	DO i = 1, Int(PDF_segments)
		PDF_V(1,i) = i * PDF_V_segment_size
	END DO

	DO i =1, (number_columns -1)
		DO j = 1, Int(PDF_segments)
			PDF_V(i+1,j) = REAL(PDF_helper(i,j))/REAL(number_rows)
		END DO
	END DO

	DEALLOCATE (PDF_helper)

	!find two max-values of the bimodal PDF function
	DO i = 1, number_columns-1
		DO j = 1, Int(PDF_upper_boundary_air)
		Max_air = max(Max_air, PDF_V(i+1,j))
		END DO

		DO k = Int(PDF_segments), Int(PDF_upper_boundary_air), -1
		Max_water =  max(Max_water, PDF_V(i+1,k))
		END DO

		DO l = 1, Int(PDF_upper_boundary_air)
			locator = locator +1
			IF (Max_air - PDF_V(i+1,l) < 0.00001) THEN
			max_air_loc = PDF_V(1,locator)
			EXIT
			END IF
		END DO

		locator = Int(PDF_segments)

		DO l = Int(PDF_segments), Int(PDF_upper_boundary_air), -1
			IF (Max_water - PDF_V(i+1,l) < 0.00001) THEN
			max_water_loc = PDF_V(1,locator)
			EXIT
			END IF
			locator = locator - 1
		END DO


	threshold(i) = max_air_loc + (max_water_loc - max_air_loc)*threshold_C
	Max_water = 0
	Max_air = 0
	locator = 0
	END DO


END SUBROUTINE Probability_V

!_______________________________________________________________________________________________________________

SUBROUTINE Basic_properties (Data_array, number_columns, number_rows, threshold,&
							 interfaces, basic_results, max_no_bubbles, frequency)

INTEGER, PARAMETER :: DBL = SELECTED_REAL_KIND(p=15)	!Precision of the Real variables
	INTEGER, INTENT (IN) :: number_columns
	INTEGER, INTENT (IN) :: number_rows
	INTEGER, INTENT (IN) :: max_no_bubbles
	INTEGER, INTENT(IN) :: frequency
	REAL(kind=DBL), INTENT (IN) :: Data_array (number_columns, number_rows)
	REAL(kind=DBL), INTENT (IN) :: threshold (number_columns-1)
	INTEGER, INTENT (INOUT) :: interfaces((number_columns-1)*2, max_no_bubbles)  
	REAL(kind=DBL), INTENT (INOUT) :: basic_results ((number_columns-1)*2, 1)
	INTEGER, ALLOCATABLE, DIMENSION  (:,:) :: instantenious_C !helper -> air=1, water =0 
	INTEGER i, j, status	
	INTEGER :: size_counter, bubble_counter, C_counter

	ALLOCATE (instantenious_C(number_columns-1, number_rows))
	instantenious_C = 0
	interfaces = 0
	basic_results = 0
	status = 0
	size_counter = 0
	bubble_counter = 0
	C_counter = 0

	DO i =1, (number_columns-1)
		DO j = 1, number_rows
			IF (Data_array((i+1),j) - threshold (i) < 0.00001) THEN
				instantenious_C(i, j) = 1
				ELSE
				instantenious_C(i, j) = 0
			END IF
		END DO
	END DO	
		

	DO i =1, (number_columns-1)					!for all probes
	C_counter =0
		IF (instantenious_C(i,1) == 0) THEN
			size_counter = 0
			bubble_counter = 0
			ELSE 
			size_counter = 1
			bubble_counter = 1
		END IF

		DO j = 1, number_rows-1
			IF (instantenious_C(i,j) == 1) THEN
			C_counter = C_counter +1
			END IF

			IF (instantenious_C(i,j) == 0 .AND. instantenious_C(i,j+1) == 1 ) THEN
				size_counter = size_counter +1
				interfaces((i*2-1), size_counter) = j+1				!change water to air
				bubble_counter = bubble_counter +1
			ELSE IF (instantenious_C(i,j) == 1 .AND. instantenious_C(i,j+1) == 0 ) THEN	
				interfaces((i*2), size_counter) = j	+1			!change air to water	
			END IF
		END DO
		interfaces((i*2), size_counter) = number_rows
		basic_results(i*2-1,1) = REAL(C_counter)/REAL(number_rows)
		basic_results(i*2, 1) = REAL(bubble_counter)*REAL(frequency)/REAL(number_rows)	
		
	END DO
	
	DEALLOCATE (instantenious_C)

END SUBROUTINE Basic_properties



!_______________________________________________________________________________________________________________

SUBROUTINE Filter (number_rows, number_columns, FFT_number, Int_spectral_analys, FFT_array, Cutoff_time_low_lower, &
				Cutoff_time_low_upper, Cutoff_time_mid_lower, Cutoff_time_mid_upper, Cutoff_time_high_lower, &
				Cutoff_time_high_upper, Filtered_array_low, Filtered_array_mid, Filtered_array_high)

INTEGER, PARAMETER :: DBL = SELECTED_REAL_KIND(p=15)	!Precision of the Real variables
INTEGER, INTENT (IN) :: number_rows
INTEGER, INTENT (IN) :: number_columns
INTEGER, INTENT (IN) :: FFT_number
INTEGER, INTENT (IN) :: Int_spectral_analys
REAL(kind=DBL), INTENT (IN) :: FFT_array (number_columns, FFT_number)
INTEGER, INTENT (IN) :: Cutoff_time_low_lower
INTEGER, INTENT (IN) :: Cutoff_time_low_upper
INTEGER, INTENT (IN) :: Cutoff_time_mid_lower
INTEGER, INTENT (IN) :: Cutoff_time_mid_upper
INTEGER, INTENT (IN) :: Cutoff_time_high_lower
INTEGER, INTENT (IN) :: Cutoff_time_high_upper
REAL(kind=DBL), INTENT (INOUT) :: Filtered_array_low (number_columns, number_rows)
REAL(kind=DBL), INTENT (INOUT) :: Filtered_array_mid (number_columns, number_rows)
REAL(kind=DBL), INTENT (INOUT) :: Filtered_array_high (number_columns, number_rows)
INTEGER :: two_times_Int_spectral_analys, four_times_Int_spectral_analys
INTEGER i, j, k, l

REAL(kind=DBL), ALLOCATABLE, DIMENSION  (:,:) :: FFT_array_helper
REAL(kind=DBL), ALLOCATABLE, DIMENSION  (:,:) :: FFT_array_inputhelper
REAL(kind=DBL), ALLOCATABLE, DIMENSION  (:) :: FFT_inputhelper
REAL(kind=DBL), ALLOCATABLE, DIMENSION  (:) :: FFT_outputhelper

ALLOCATE (FFT_array_helper (number_columns, 2*FFT_number))
ALLOCATE (FFT_array_inputhelper (number_columns, 2*FFT_number))
ALLOCATE (FFT_inputhelper (2*FFT_number))
ALLOCATE (FFT_outputhelper (2*FFT_number))
FFT_inputhelper = 0
FFT_array_helper = 0
FFT_array_inputhelper = 0
FFT_inputhelper = 0

	Do k = 2, number_columns
		Do j = 1, FFT_number
			FFT_array_helper (k,j) = FFT_array (k,j)
		END DO
	END DO

	Do k = 2, number_columns
		Do j = FFT_number + 1, (2*FFT_number)
			FFT_array_helper (k,j) = 0 
		END DO
	END DO


two_times_Int_spectral_analys = Int_spectral_analys + Int_spectral_analys
four_times_Int_spectral_analys = two_times_Int_spectral_analys + two_times_Int_spectral_analys


	! lower frequency (mean)
	Do k = 2, number_columns
		DO J=1,two_times_Int_spectral_analys
		FFT_array_inputhelper(k, 2*J-1)= FFT_array_helper(k, J)
		FFT_array_inputhelper(k, 2*J)= FFT_array_helper(k, J + two_times_Int_spectral_analys)
		END DO
	END DO


Do k = 2, number_columns
	Do j = 1, (2*FFT_number)
	FFT_inputhelper (j) = FFT_array_inputhelper(k, j)
	END DO

	CALL FOUR1 (FFT_inputhelper, two_times_Int_spectral_analys, 1)


	DO l=1,Int_spectral_analys
	  IF (l .EQ. 1.AND.Cutoff_time_low_lower > 0) THEN
	    FFT_inputhelper(1)=0.
	    FFT_inputhelper(2)=0.
	  ELSE IF (l < Cutoff_time_low_lower) THEN
	    FFT_inputhelper(2*l-1) = 0.
	    FFT_inputhelper(2*l) = 0.
	    FFT_inputhelper(four_times_Int_spectral_analys -2*l +3) =0.
	    FFT_inputhelper(four_times_Int_spectral_analys -2*l +4) =0.
	  ENDIF
	  IF (l > Cutoff_time_low_upper) THEN
	    FFT_inputhelper(2*l-1)=0.
	    FFT_inputhelper(2*l)=0.
	    FFT_inputhelper(four_times_Int_spectral_analys-2*l+3)=0.
	    FFT_inputhelper(four_times_Int_spectral_analys-2*l+4)=0.
	  ENDIF
	  IF (Cutoff_time_low_upper > Int_spectral_analys) THEN
	    FFT_inputhelper(two_times_Int_spectral_analys+1)=0.
	    FFT_inputhelper(two_times_Int_spectral_analys+2)=0.
	  ENDIF

	ENDDO
	
	CALL FOUR1(FFT_inputhelper,two_times_Int_spectral_analys, -1)


	DO l=1,two_times_Int_spectral_analys
	 FFT_outputhelper(l) = FFT_inputhelper (2*l-1)/REAL(two_times_Int_spectral_analys)
	 FFT_outputhelper(l + two_times_Int_spectral_analys) = FFT_inputhelper (2*l)/REAL(two_times_Int_spectral_analys)
	END DO

	Do j = 1, (2*FFT_number)
	FFT_array_inputhelper(k, j) = FFT_outputhelper (j)
	END DO
	
END DO


DO k = 2, number_columns
	Do j = 1, number_rows
		Filtered_array_low (k, j) = FFT_array_inputhelper(k, j)

END DO
	END DO

	FFT_array_inputhelper = 0
	FFT_inputhelper = 0


! mid (slow fluctuations)

	Do k = 2, number_columns
		DO J=1,two_times_Int_spectral_analys
		FFT_array_inputhelper(k, 2*J-1)= FFT_array_helper(k, J)
		FFT_array_inputhelper(k, 2*J)= FFT_array_helper(k, J + two_times_Int_spectral_analys)
		END DO
	END DO

Do k = 2, number_columns
	Do j = 1, (2*FFT_number)
	FFT_inputhelper (j) = FFT_array_inputhelper(k, j)
	END DO

	CALL FOUR1 (FFT_inputhelper, two_times_Int_spectral_analys, 1)


	DO l=1,Int_spectral_analys
	  IF (l .EQ. 1.AND.Cutoff_time_mid_lower > 0) THEN
	    FFT_inputhelper(1)=0.
	    FFT_inputhelper(2)=0.
	  ELSE IF (l < Cutoff_time_mid_lower) THEN
	    FFT_inputhelper(2*l-1) = 0.
	    FFT_inputhelper(2*l) = 0.
	    FFT_inputhelper(four_times_Int_spectral_analys -2*l +3) =0.
	    FFT_inputhelper(four_times_Int_spectral_analys -2*l +4) =0.
	  ENDIF
	  IF (l > Cutoff_time_mid_upper) THEN
	    FFT_inputhelper(2*l-1)=0.
	    FFT_inputhelper(2*l)=0.
	    FFT_inputhelper(four_times_Int_spectral_analys-2*l+3)=0.
	    FFT_inputhelper(four_times_Int_spectral_analys-2*l+4)=0.
	  ENDIF
	  IF (Cutoff_time_mid_upper > Int_spectral_analys) THEN
	    FFT_inputhelper(two_times_Int_spectral_analys+1)=0.
	    FFT_inputhelper(two_times_Int_spectral_analys+2)=0.
	  ENDIF

	ENDDO
	
	CALL FOUR1(FFT_inputhelper,two_times_Int_spectral_analys, -1)


	DO l=1,two_times_Int_spectral_analys
	 FFT_outputhelper(l) = FFT_inputhelper (2*l-1)/REAL(two_times_Int_spectral_analys)
	 FFT_outputhelper(l + two_times_Int_spectral_analys) = FFT_inputhelper (2*l)/REAL(two_times_Int_spectral_analys)
	END DO

	Do j = 1, (2*FFT_number)
	FFT_array_inputhelper(k, j) = FFT_outputhelper (j)
	END DO
	
END DO


DO k = 2, number_columns
	Do j = 1, number_rows
		Filtered_array_mid (k, j) = FFT_array_inputhelper(k, j)

END DO
	END DO

	FFT_array_inputhelper = 0
	FFT_inputhelper = 0



! high (high fluctuations)

	Do k = 2, number_columns
		DO J=1,two_times_Int_spectral_analys
		FFT_array_inputhelper(k, 2*J-1)= FFT_array_helper(k, J)
		FFT_array_inputhelper(k, 2*J)= FFT_array_helper(k, J + two_times_Int_spectral_analys)
		END DO
	END DO

Do k = 2, number_columns
	Do j = 1, (2*FFT_number)
	FFT_inputhelper (j) = FFT_array_inputhelper(k, j)
	END DO

	CALL FOUR1 (FFT_inputhelper, two_times_Int_spectral_analys, 1)


	DO l=1,Int_spectral_analys
	  IF (l .EQ. 1.AND.Cutoff_time_high_lower > 0) THEN
	    FFT_inputhelper(1)=0.
	    FFT_inputhelper(2)=0.
	  ELSE IF (l < Cutoff_time_high_lower) THEN
	    FFT_inputhelper(2*l-1) = 0.
	    FFT_inputhelper(2*l) = 0.
	    FFT_inputhelper(four_times_Int_spectral_analys -2*l +3) =0.
	    FFT_inputhelper(four_times_Int_spectral_analys -2*l +4) =0.
	  ENDIF
	  IF (l > Cutoff_time_high_upper) THEN
	    FFT_inputhelper(2*l-1)=0.
	    FFT_inputhelper(2*l)=0.
	    FFT_inputhelper(four_times_Int_spectral_analys-2*l+3)=0.
	    FFT_inputhelper(four_times_Int_spectral_analys-2*l+4)=0.
	  ENDIF
	  IF (Cutoff_time_high_upper > Int_spectral_analys) THEN
	    FFT_inputhelper(two_times_Int_spectral_analys+1)=0.
	    FFT_inputhelper(two_times_Int_spectral_analys+2)=0.
	  ENDIF

	ENDDO
	
	CALL FOUR1(FFT_inputhelper,two_times_Int_spectral_analys, -1)


	DO l=1,two_times_Int_spectral_analys
	 FFT_outputhelper(l) = FFT_inputhelper (2*l-1)/REAL(two_times_Int_spectral_analys)
	 FFT_outputhelper(l + two_times_Int_spectral_analys) = FFT_inputhelper (2*l)/REAL(two_times_Int_spectral_analys)
	END DO

	Do j = 1, (2*FFT_number)
	FFT_array_inputhelper(k, j) = FFT_outputhelper (j)
	END DO
END DO


DO k = 2, number_columns
	Do j = 1, number_rows
		Filtered_array_high (k, j) = FFT_array_inputhelper(k, j)

END DO
	END DO

DEALLOCATE (FFT_array_helper)
DEALLOCATE (FFT_array_inputhelper)
DEALLOCATE (FFT_inputhelper)
DEALLOCATE (FFT_outputhelper)

		WRITE (*,*) 'Finished filtering for triple-decomposition of raw-signal! '

END SUBROUTINE Filter



!_______________________________________________________________________________________________________________
		
SUBROUTINE FOUR1(D,NN, ISIGN)

!	based on four1(sub-routine) Numerical recipes in FORTRAN by Press 
!       et.al. 1992
INTEGER, PARAMETER :: DBL = SELECTED_REAL_KIND(p=15)	!Precision of the Real variables
	INTEGER ISIGN,NN
	REAL(kind=DBL) D(2*NN)
	INTEGER I,ISTEP,J,M,MMAX,N
	REAL(kind=DBL) TEMPI,TEMPR
	REAL(kind=DBL) THETA,WI,WPI,WPR,WR,WTEMP

	N=2*NN
	J=1

	DO I=1,N,2
	  IF (J.GT.I) THEN
	    TEMPR=D(J)
	    TEMPI=D(J+1)
	    D(J)=D(I)
            D(J+1)=D(I+1)
	    D(I)=TEMPR
	    D(I+1)=TEMPI
	  ENDIF
	  M=N/2
1	  IF ((M.GE.2).AND.(J.GT.M)) THEN
	    J=J-M
	    M=M/2
	    GOTO 1
	  ENDIF
	  J=J+M
	ENDDO

	MMAX=2
2	IF (N.GT.MMAX) THEN
	  ISTEP=2*MMAX
	  THETA=6.28318530717959D0/(ISIGN*MMAX)
	  WPR=-2.D0*DSIN(0.5D0*THETA)**2
	  WPI=DSIN(THETA)
	  WR=1.D0
	  WI=0.D0
	  DO M=1,MMAX,2
	    DO I=M,N,ISTEP
	      J=I+MMAX
	      TEMPR=SNGL(WR)*D(J)-SNGL(WI)*D(J+1)
	      TEMPI=SNGL(WR)*D(J+1)+SNGL(WI)*D(J)
	      D(J)=D(I)-TEMPR
	      D(J+1)=D(I+1)-TEMPI
	      D(I)=D(I)+TEMPR
	      D(I+1)=D(I+1)+TEMPI
	    ENDDO 
	    WTEMP=WR
	    WR=WR*WPR-WI*WPI+WR
	    WI=WI*WPR+WTEMP*WPI+WI
	  ENDDO 
	  MMAX=ISTEP
	  GOTO 2
	ENDIF

	RETURN
	END SUBROUTINE FOUR1


