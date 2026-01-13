BASE   = ../..
NAME   = $(shell basename `pwd`)
OUTPUT = lib$(NAME).so

OPTIM = debug

DEPEND = Core Algo extern_CMA extern_liblbfgs

CERES = 0
NLOPT = 0
IPOPT = 0

SRCS = $(shell find . -maxdepth 1 -name '*.cpp' )
OBJS = $(SRCS:%.cpp=%.o)

include $(BASE)/_make/generic.mk
