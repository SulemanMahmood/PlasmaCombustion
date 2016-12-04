CHARMC = /home/kchndrs2/charm/netlrts-linux-x86_64/bin/charmc
OBJS = run.o

all: run

run: $(OBJS)
	$(CHARMC) -O3 -language -charm++ -o run $(OBJS)

run.decl.h: run.ci
	$(CHARMC) run.ci

clean:
	rm -f *.decl.h *.def.h conv-host *.o run charmrun

cleanp:
	rm -f *.sts *.gz *.projc *.topo *.out

run.o: run.C main.decl.h
	$(CHARMC) -c run.C

test: all
	./charmrun run +p32
