def:
	@test -d obj || mkdir obj
	@cp src/Makefile obj
	@cd obj && $(MAKE) COMPENV="-DM2MSTD -DM2LSTD -DL2LSTD -DKSTD" EXE=../fmm.x

all:
	$(MAKE) -C obj clean EXE=../fmm.x-rku
	@test -d obj || mkdir obj
	@cp src/Makefile obj
	@cd obj && $(MAKE) COMPENV="-DM2MROT -DM2LROT -DL2LROT -DKSTD" EXE=../fmm.x-rku

	$(MAKE) -C obj clean EXE=../fmm.x-rmu
	@test -d obj || mkdir obj
	@cp src/Makefile obj
	@cd obj && $(MAKE) COMPENV="-DM2MROT -DM2LROT -DL2LROT -DMSTD" EXE=../fmm.x-rmu

	$(MAKE) -C obj clean EXE=../fmm.x-rkk
	@test -d obj || mkdir obj
	@cp src/Makefile obj
	@cd obj && $(MAKE) COMPENV="-DM2MROT -DM2LROT -DL2LROT -DKSCA" EXE=../fmm.x-rkk

	$(MAKE) -C obj clean EXE=../fmm.x-rmm
	@test -d obj || mkdir obj
	@cp src/Makefile obj
	@cd obj && $(MAKE) COMPENV="-DM2MROT -DM2LROT -DL2LROT -DMSCA" EXE=../fmm.x-rmm

	$(MAKE) -C obj clean EXE=../fmm.x-rmk
	@test -d obj || mkdir obj
	@cp src/Makefile obj
	@cd obj && $(MAKE) COMPENV="-DM2MROT -DM2LROT -DL2LROT -DMKSCA" EXE=../fmm.x-rmk

	$(MAKE) -C obj clean EXE=../fmm.x-std
	@test -d obj || mkdir obj
	@cp src/Makefile obj
	@cd obj && $(MAKE) COMPENV="-DM2MSTD -DM2LSTD -DL2LSTD -DKSTD" EXE=../fmm.x-std

rku:
	@test -d obj || mkdir obj
	@cp src/Makefile obj
	@cd obj && $(MAKE) COMPENV="-DM2MROT -DM2LROT -DL2LROT -DKSTD" EXE=../fmm.x-rku
rmk:
	@test -d obj || mkdir obj
	@cp src/Makefile obj
	@cd obj && $(MAKE) COMPENV="-DM2MROT -DM2LROT -DL2LROT -DMKSCA" EXE=../fmm.x-rmk
rkk:
	@test -d obj || mkdir obj
	@cp src/Makefile obj
	@cd obj && $(MAKE) COMPENV="-DM2MROT -DM2LROT -DL2LROT -DKSCA" EXE=../fmm.x-rkk
rrs:
	@test -d obj || mkdir obj
	@cp src/Makefile obj
	@cd obj && $(MAKE) COMPENV="-DM2MROT -DM2LSTD -DL2LROT -DKSCA" EXE=../fmm.x-rrs
rmm:
	@test -d obj || mkdir obj
	@cp src/Makefile obj
	@cd obj && $(MAKE) COMPENV="-DM2MROT -DM2LROT -DL2LROT -DMSCA" EXE=../fmm.x-rmm
std:
	@test -d obj || mkdir obj
	@cp src/Makefile obj
	@cd obj && $(MAKE) COMPENV="-DM2MSTD -DM2LSTD -DL2LSTD -DKSTD" EXE=../fmm.x-std

clean:
	$(MAKE) -C obj clean
