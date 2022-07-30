 FC = ifort 
 TARGET = testFMM.out

 OBJDIR = ./

 SRCS = testFMM.f90
 OBJS = testFMM.o

 #FFLAGS = -check all
 #FFLAGS = -O2 -xHost -ipo

.SUFFIXES: .f90

$(TARGET): $(OBJS)
	$(FC) -o $@ $^ $(FFLAGS) 
	rm -f $(OBJS) ./*.mod ./*.o

%.o %.mod: %.f90
	$(FC) -o $@ -c $<

all: $(TARGET)

clean:
	rm -f $(OBJS) ./*.mod
