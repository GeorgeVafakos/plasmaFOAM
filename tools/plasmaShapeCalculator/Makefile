# Makefile

CC = gcc
GFLAGS = -Wall -lm

SOURCES = plasmaShapeCalculator.c
OBJECTS = $(SOURCES:.c=.o)
EXECUTABLE = plasmaShapeCalculator

PROGRESS = \033[1;91mProgress:\033[0m

PERCENT = 0.01

execname: $(SOURCES)
	@$(CC) $(SOURCES) $(GFLAGS) -o $(EXECUTABLE)
	

	@echo "\033[1mStart compiling"
	@echo -ne '$(PROGRESS) \033[1m  ..........                                                                                              \033[91m[10%]\033[0m\r'
	@sleep $(PERCENT)
	@echo -ne '$(PROGRESS) \033[1m  ...............                                                                                         \033[91m[10%]\033[0m\r'
	@sleep $(PERCENT)
	@echo -ne '$(PROGRESS) \033[1m  ....................                                                                                    \033[91m[20%]\033[0m\r'
	@sleep $(PERCENT)
	@echo -ne '$(PROGRESS) \033[1m  .........................                                                                               \033[91m[20%]\033[0m\r'
	@sleep $(PERCENT)
	@echo -ne '$(PROGRESS) \033[1m  ..............................                                                                          \033[91m[30%]\033[0m\r'
	@sleep $(PERCENT)
	@echo -ne '$(PROGRESS) \033[1m  ...................................                                                                     \033[91m[30%]\033[0m\r'
	@sleep $(PERCENT)
	@echo -ne '$(PROGRESS) \033[1m  ........................................                                                                \033[91m[40%]\033[0m\r'
	@sleep $(PERCENT)
	@echo -ne '$(PROGRESS) \033[1m  .............................................                                                           \033[91m[40%]\033[0m\r'
	@sleep $(PERCENT)
	@echo -ne '$(PROGRESS) \033[1m  ..................................................                                                      \033[91m[50%]\033[0m\r'
	@sleep $(PERCENT)
	@echo -ne '$(PROGRESS) \033[1m  .......................................................                                                 \033[91m[50%]\033[0m\r'
	@sleep $(PERCENT)
	@echo -ne '$(PROGRESS) \033[1m  ............................................................                                            \033[91m[60%]\033[0m\r'
	@sleep $(PERCENT)
	@echo -ne '$(PROGRESS) \033[1m  .................................................................                                       \033[91m[60%]\033[0m\r'
	@sleep $(PERCENT)
	@echo -ne '$(PROGRESS) \033[1m  ......................................................................                                  \033[91m[70%]\033[0m\r'
	@sleep $(PERCENT)
	@echo -ne '$(PROGRESS) \033[1m  ...........................................................................                             \033[91m[70%]\033[0m\r'
	@sleep $(PERCENT)
	@echo -ne '$(PROGRESS) \033[1m  ................................................................................                        \033[91m[80%]\033[0m\r'
	@sleep $(PERCENT)
	@echo -ne '$(PROGRESS) \033[1m  .....................................................................................                   \033[91m[80%]\033[0m\r'
	@sleep $(PERCENT)
	@echo -ne '$(PROGRESS) \033[1m  ..........................................................................................              \033[91m[90%]\033[0m\r'
	@sleep $(PERCENT)
	@echo -ne '$(PROGRESS) \033[1m  ...............................................................................................         \033[91m[90%]\033[0m\r'
	@sleep $(PERCENT)
	@echo -ne '$(PROGRESS) \033[1m  ....................................................................................................    \033[91m[100%]\033[0m\r'
	@echo -ne '\n'
	@echo "\033[1mExecutable created"
	@echo "Compilation was successful!"

clean:
	rm -rf $(EXECUTABLE)
	@echo "\033[1mFiles removed!"
