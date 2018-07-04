OBJ=ampt.o amptsub.o art1f.o \
hijing1.383_ampt.o hipyset1.35.o zpc.o linana.o
F77=gfortran -ffree-line-length-0

all: ampt clean

ampt: $(OBJ)
	$(F77) -o ampt -O $(OBJ)

ampt.o: ampt.f90
	$(F77)  -c -O ampt.f90

amptsub.o: amptsub.f90
	$(F77)  -c -O amptsub.f90

art1f.o: art1f.f90
	$(F77)  -c -O art1f.f90

hijing1.383_ampt.o: hijing1.383_ampt.f90
	$(F77)  -c -O hijing1.383_ampt.f90

hipyset1.35.o: hipyset1.35.f90
	$(F77)  -c -O hipyset1.35.f90

zpc.o: zpc.f90
	$(F77)  -c -O zpc.f90

linana.o: linana.f90
	$(F77)  -c -O linana.f90

.PHONY:clean
clean:
	$(RM) $(OBJ)

.PHONY:distclean
distclean:
	$(RM) -r ampt
