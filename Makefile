# ######################################################################
# Makefile for seismic-gf-uncertainty project
# 
# Purpose: Compiles Fortran90 modules and example program for evaluation 
#          of approximate covariance matrices (ACF, AXCF, SACF, SAXCF).
# 
# Author:  Miroslav HALLO
# Date:    2026/04
# License: GNU GPL v3
# ######################################################################
#
# Usage:
#    > make            # Compile all binary files
#    > make clean      # Remove all binary files and modules
#    > make run        # Compile all and run the main program
#    > make clean run  # Remove all, Compile all, and run the program
#
# ######################################################################

# Compiler
F90 = gfortran

# Flags
FFLAGS = -O2 -fcheck=all

# Program name
TARGET = run_example

# Path to codes
SRCDIR = src

# Oject list
OBJ = $(SRCDIR)/nr.o $(SRCDIR)/approxc.o $(SRCDIR)/example.o

# ----------------------------------------------------------------------

# Compile main program
$(TARGET): $(OBJ)
	$(F90) $(FFLAGS) -o $@ $(OBJ)

# Object for Fortran 77 routine from Numerical Recipes
$(SRCDIR)/nr.o: $(SRCDIR)/nr.for
	$(F90) $(FFLAGS) -c $< -o $@

# Object for modul with ACF, AXCF, SACF, SAXCF
$(SRCDIR)/approxc.o: $(SRCDIR)/approxc.f90
	$(F90) $(FFLAGS) -c $< -o $@ -J$(SRCDIR)

# Object for main program
$(SRCDIR)/example.o: $(SRCDIR)/example.f90 $(SRCDIR)/approxc.o
	$(F90) $(FFLAGS) -c $< -o $@ -I$(SRCDIR)

# ----------------------------------------------------------------------
# Clean after compilation
clean:
	rm -f $(SRCDIR)/*.o $(SRCDIR)/*.mod $(TARGET)

# ----------------------------------------------------------------------
# Run the main program
run: $(TARGET)
	./$(TARGET)


.PHONY: clean run