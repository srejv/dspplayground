INCLUDES = -I./include -I/usr/local/include
LIBS	 = -L./lib -L/usr/local/lib -lportaudio -lm
CC = gcc

# säkert inte såhär man gör 
all:
	make helloring
	make hellotable
	make hellodelay

helloring:	
	$(CC) -o bin/helloring src/helloring.c $(INCLUDES) $(LIBS)

hellotable:	
	$(CC) -o bin/hellotable src/hellotable.c src/wave.c src/gtable.c $(INCLUDES) $(LIBS)

hellodelay:
	$(CC) -o bin/hellodelay src/hellodelay.c src/delay.c src/gtable.c src/wave.c $(INCLUDES) $(LIBS)

helloflanger:
	$(CC) -o bin/helloflanger src/helloflanger.c src/bk.c $(INCLUDES) $(LIBS)

hellopitch:
	$(CC) -o bin/hellopitch src/hellopitch.c src/bk.c $(INCLUDES) $(LIBS)

helloharmonizer:
	$(CC) -o bin/helloharmonizer src/helloharmonizer.c src/bk.c $(INCLUDES) $(LIBS)

helloenv:
	$(CC) -o bin/envtest src/envtest.c $(INCLUDES) $(LIBS)

hellomtap:
	$(CC) -o bin/hellomtap src/hellomtap.c src/mtapdelay.c $(INCLUDES) $(LIBS)

hellodist:
	$(CC) -o bin/hellodist src/hellodist.c $(INCLUDES) $(LIBS)