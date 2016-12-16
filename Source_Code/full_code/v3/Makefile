CHARMC = charmc

all:
	$(CHARMC) run.ci;
	$(CHARMC) -c run.C;
	$(CHARMC) -o run run.o -tracemode projections

clean:
	rm Output.csv run.o run charmrun *.sum *.out *.projrc *.sts run.decl.h run.def.h
	
test: all
	./charmrun run +p4 ++local