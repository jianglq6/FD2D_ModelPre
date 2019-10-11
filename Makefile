#=======================================
#
#  For Lebedev grid
#
#  2019.09
#
#==========================================

#=======================================================
# Directory
#=======================================================
SRCDIR := ./src
BINDIR := ./bin
OBJDIR := ./OBJ
skeldirs := $(OBJDIR) $(BINDIR)

#=======================================================
# Compiler and Option
#=======================================================
CC := gcc
CFLAGS := -c -g -I $(SRCDIR) #-std=c99
#CFLAGS := -c

#=======================================================
# Source File Names
#=======================================================

SRC := Elastic2d_tti.c elastic2d_src.c elastic2d_lebedev.c staggered_fd_coef.c elastic2d_filter.c elastic2d_math.c elastic2d_abs_exp.c  elastic2d_stf.c read_config_para.c share_param.c write_snapshots.c write_seismograms.c 
INC := elastic2d_src.h elastic2d_lebedev.h  staggered_fd_coef.h elastic2d_filter.h elastic2d_math.h elastic2d_abs_exp.h  elastic2d_stf.h read_config_para.h share_param.h write_snapshots.h Elastic2d_tti.h write_seismograms.h
#OBJ := $(SRC:%.c=%.o)
OBJ := $(foreach file,$(SRC),$(OBJDIR)/$(file:.c=.o))

EXE := Lebedev_tti
#=======================================================
#Searching Path
#=======================================================
vpath %.c $(SRCDIR)
vpath %.h $(SRCDIR)

#=======================================================
# Target
#=======================================================
PHONYLIST = skel all slover
.PHONY : $(PHINYLIST)

all: skel solver

solver: $(BINDIR)/$(EXE)

skel:
	@echo $(OBJ)
	@mkdir -p $(skeldirs)

$(BINDIR)/$(EXE) : $(OBJ)
	$(CC) -o $@ $(^) -lm
	@echo "----- DONE -----"
	@echo

$(OBJDIR)/%.o : %.c
	$(CC) $(CFLAGS) $< -o $(@) -lm

.PHONY : clean

RM := rm
clean:
	@echo "cleaning..."
	$(RM) -r $(BINDIR)/*
	$(RM) -R $(OBJDIR)/*
cleanexe:
	$(RM) -f $(BINDIR)/*
cleanall: cleanexe clean
distclean: cleanall
	$(RM) -r $(BINDIR)/*
	$(RM) -R $(OBJDIR)/*
