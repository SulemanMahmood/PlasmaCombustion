CHARMC = charmc

all:
	$(CHARMC) run.ci;
	$(CHARMC) -c run.C;
	$(CHARMC) -o run run.o
	
test: all
	./charmrun run +p2 ++local