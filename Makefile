#=======================================================
#
#  For Lebedev grid
#
#  2019.09
#
#=======================================================

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

SRC := Elastic2d_FD.c elastic2d_src.c elastic2d_lebedev.c elastic2d_staggered.c staggered_fd_coef.c elastic2d_filter.c elastic2d_math.c elastic2d_abs_exp.c  elastic2d_stf.c read_config_para.c share_param.c write_snapshots.c write_seismograms.c pre_model_prepare.c read_interface_material.c pre_para_assign.c pre_cal_dip_area.c pre_media_parameterization.c
INC := elastic2d_src.h elastic2d_lebedev.h  elastic2d_staggered.h staggered_fd_coef.h elastic2d_filter.h elastic2d_math.h elastic2d_abs_exp.h  elastic2d_stf.h read_config_para.h share_param.h write_snapshots.h write_seismograms.h pre_interface_struct.h read_interface_material.h pre_para_assign.h pre_cal_dip_area.h pre_media_parameterization.h pre_model_prepare.h
#OBJ := $(SRC:%.c=%.o)
OBJ := $(foreach file,$(SRC),$(OBJDIR)/$(file:.c=.o))

EXE := FD2D_modeling
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
