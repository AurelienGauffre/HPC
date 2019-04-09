# Options complementaires pour l'edition de liens
LDFLAGS=

# Options de compilation
CXXFLAGS=-Wall

# Choix du compilateur
CXX=g++

# Executables a creer
SRC = $(wildcard *.cxx)
OBJ = $(SRC:.cxx=.o)
EXEC = string

#
FLAGS = -DVERBOSE=0
all : $(EXEC)

# Regle de creation de binaire
string: string_equation.o
	$(CXX) $(CXXFLAGS) -o $@ $+

# Regle de creation d'objets
%.o: %.cxx
	$(CXX) $(FLAGS) $(CXXFLAGS) -c $<

# Nettoye les fichiers de compilations intermediaires
clean:
	rm -f *.o *~

# Nettoye l'integralite des fichiers compiles
mrproper: clean
	rm $(EXEC)
