CHARMC = /home/kchndrs2/charm/netlrts-linux-x86_64/bin/charmc

all:
	$(CHARMC) run.ci;
	$(CHARMC) -c run.C;
	$(CHARMC) -o run run.o -tracemode projections

clean:
	rm run.o run charmrun *.sum *.out *.projrc *.sts run.decl.h run.def.h

test: all
	./charmrun run +p10
