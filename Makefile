F77 = gfortran

SRC_PATH   := ./src
SRCS       += $(wildcard $(SRC_PATH)/*.f)

all: clean build

.PHONY:build
build:
	@$(F77) $(SRCS) -o ampt

.PHONY:clean
clean:
	@$(RM) ampt
