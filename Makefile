# --- Konfiguracja kompilatora ---

# Kompilator C++
CXX = g++

# Flagi kompilatora:
# -Wall = włącz wszystkie ostrzeżenia
# -g = dołącz symbole debugowania (dla gdb)
# -std=c++11 = użyj standardu C++11 (wymagane m.in. dla NAN, chrono, random)
# -MMD -MP = flagi do automatycznego generowania zależności od plików .h
CXXFLAGS = -Wall -g -std=c++11 -MMD -MP

# Flagi dla linkera:
# -lm = dołącz bibliotekę matematyczną (dla <cmath>)
LDFLAGS = -lm

# --- Konfiguracja projektu ---

# Nazwa docelowego pliku wykonywalnego
TARGET = optimizer

# Znajdź wszystkie pliki .cpp w bieżącym katalogu
SOURCES = $(wildcard *.cpp)

# Wygeneruj listę plików obiektowych (.o) na podstawie plików źródłowych (.cpp)
OBJECTS = $(SOURCES:.cpp=.o)

# Wygeneruj listę plików zależności (.d)
DEPS = $(SOURCES:.cpp=.d)

# --- Reguły ---

# Reguła domyślna (wywoływana przez 'make')
all: $(TARGET)

# Reguła do linkowania programu
$(TARGET): $(OBJECTS)
	@echo "Łączenie programu..."
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJECTS) $(LDFLAGS)
	@echo "Budowanie zakończone: $(TARGET)"

# Reguła wzorcowa do kompilacji plików .cpp do .o
# $< to nazwa pierwszej zależności (pliku .cpp)
# $@ to nazwa celu (pliku .o)
%.o: %.cpp
	@echo "Kompilowanie $<..."
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Reguła do czyszczenia skompilowanych plików
clean:
	@echo "Czyszczenie plików..."
	rm -f $(TARGET) $(OBJECTS) $(DEPS)
	@echo "Gotowe."

# Reguła do uruchomienia programu
run: all
	./$(TARGET)

# Dołącza wygenerowane pliki zależności.
# Znak '-' na początku ignoruje błąd, jeśli plik .d jeszcze nie istnieje.
-include $(DEPS)

# Definiuje cele, które nie są plikami
.PHONY: all clean run