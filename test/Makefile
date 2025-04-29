



phantom_dats = phantom1.dat \
           phantom2.dat \
           phantom3.dat

phantom_refs = phantom1_ref.dat \
              phantom2_ref.dat \
              phantom3_ref.dat

phantom_tests = phantom1_test \
              phantom2_test \
              phantom3_test

threshhold = 0.0
all: test


refs: $(phantom_refs) 
test: $(phantom_tests)
phantoms: $(phantom_dats)
	
phantom1.dat: phantom1.pro
	@echo 'phantom1' | idl
phantom1_ref.dat: linuxtest.pro phantom1.pro phantom1.dat
	@echo 'linuxtest,"phantom1.dat",/saveref,referencefile="phantom1_ref.dat"' | idl
phantom1_test: 
	@echo 'linuxtest,"phantom1.dat",mse,referencefile="phantom1_ref.dat" & print & print,"Phantom1 MSE: ",mse & if mse gt $(threshhold) then exit,status=-1 & print,"Passed"&print' | idl

phantom2.dat: phantom2.pro
	@echo 'phantom2' | idl
phantom2_ref.dat: linuxtest.pro phantom2.pro phantom2.dat
	@echo 'linuxtest,"phantom2.dat",/saveref,referencefile="phantom2_ref.dat"' | idl
phantom2_test: 
	@echo 'linuxtest,"phantom2.dat",mse,referencefile="phantom2_ref.dat" & print & print,"Phantom2 MSE: ",mse & if mse gt $(threshhold) then exit,status=-2 & print,"Passed"&print' | idl

phantom3.dat: phantom3.pro
	@echo 'phantom3' | idl
phantom3_ref.dat: linuxtest.pro phantom3.pro phantom3.dat
	@echo 'linuxtest,"phantom3.dat",/saveref,referencefile="phantom3_ref.dat"' | idl
phantom3_test: 
	@echo 'linuxtest,"phantom3.dat",mse,referencefile="phantom3_ref.dat" & print & print,"Phantom3 MSE: ",mse & if mse gt $(threshhold) then exit,status=-3 & print,"Passed"&print' | idl
	
